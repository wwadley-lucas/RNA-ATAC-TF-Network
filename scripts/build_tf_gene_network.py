#!/usr/bin/env python3
"""
TF-Gene Regulatory Network Builder
==================================
Integrates TOBIAS TF binding, ATAC accessibility, and RNA expression
to build TF→gene regulatory edges.

Part of RNA-ATAC-TF-Network pipeline.
https://github.com/wwadley-lucas/RNA-ATAC-TF-Network

Usage:
    python build_tf_gene_network.py --config config.yaml --contrast NvT
    python build_tf_gene_network.py --config config.yaml --contrast CvNC
"""

__version__ = "1.0.0"

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats


def load_config(config_path: str) -> dict:
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def resolve_path(path: str, base_dir: Path) -> Path:
    """Resolve relative paths from config file location."""
    p = Path(path)
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()


def load_bindetect_results(path: Path) -> pd.DataFrame:
    """Load TOBIAS BINDetect results with TF-level differential binding."""
    print(f"  Loading BINDetect results from {path.name}...")
    df = pd.read_csv(path, sep='\t')
    # Clean up column names (remove quotes)
    df.columns = [c.strip('"') for c in df.columns]
    print(f"    Loaded {len(df)} TFs with {len(df.columns)} columns")
    return df


def load_tf_binding_matrix(path: Path) -> pd.DataFrame:
    """Load TF binding count matrix (TF x samples)."""
    print(f"  Loading TF binding matrix from {path.name}...")
    df = pd.read_csv(path, index_col=0)
    # Clean TF names (remove C_ prefix if present)
    df.index = [tf.replace('C_', '') if tf.startswith('C_') else tf for tf in df.index]
    print(f"    Loaded {len(df)} TFs x {len(df.columns)} samples")
    return df


def load_peak_annotations(path: Path) -> pd.DataFrame:
    """Load peak annotations with gene mapping."""
    print(f"  Loading peak annotations from {path.name}...")
    df = pd.read_csv(path, sep='\t')
    # Create peak_id from coordinates
    df['peak_id'] = df['seqnames'].astype(str) + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    print(f"    Loaded {len(df)} peak annotations")
    return df


def load_atac_gene_counts(path: Path) -> pd.DataFrame:
    """Load gene-level ATAC accessibility counts."""
    print(f"  Loading ATAC gene counts from {path.name}...")
    df = pd.read_csv(path, sep='\t', index_col=0)
    print(f"    Loaded {len(df)} genes x {len(df.columns)} samples")
    return df


def load_rna_counts(path: Path) -> pd.DataFrame:
    """Load RNA expression count matrix."""
    print(f"  Loading RNA counts from {path.name}...")
    df = pd.read_csv(path, index_col=0)
    # Drop Description column if present
    if 'Description' in df.columns:
        df = df.drop(columns=['Description'])
    print(f"    Loaded {len(df)} genes x {len(df.columns)} samples")
    return df


def load_metadata(path: Path) -> pd.DataFrame:
    """Load sample metadata."""
    print(f"  Loading metadata from {path.name}...")
    df = pd.read_csv(path)
    print(f"    Loaded {len(df)} samples")
    return df


def create_peak2gene_mapping(peak_annot: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Create peak to gene mapping with distance weights.

    Returns DataFrame with columns: peak_id, gene, dist_to_tss, weight, annotation_type
    """
    print("\nCreating peak-to-gene mapping...")
    p2g_cfg = cfg['peak2gene']

    # Extract relevant columns
    mapping = peak_annot[['peak_id', 'geneId', 'distanceToTSS', 'annotation']].copy()
    mapping.columns = ['peak_id', 'gene', 'dist_to_tss', 'annotation']

    # Filter by annotation type
    keep_mask = pd.Series([True] * len(mapping))

    if p2g_cfg.get('use_promoter_only', False):
        keep_mask = mapping['annotation'].str.contains('Promoter', case=False, na=False)
    else:
        # Build mask based on config
        annotation_masks = []

        # Always include promoters
        annotation_masks.append(mapping['annotation'].str.contains('Promoter', case=False, na=False))

        if p2g_cfg.get('include_intron', True):
            annotation_masks.append(mapping['annotation'].str.contains('Intron', case=False, na=False))

        if p2g_cfg.get('include_distal', True):
            annotation_masks.append(mapping['annotation'].str.contains('Distal|Intergenic', case=False, na=False))

        keep_mask = pd.concat(annotation_masks, axis=1).any(axis=1)

    mapping = mapping[keep_mask].copy()

    # Filter by max distance
    max_dist = p2g_cfg.get('max_distance', 100000)
    mapping = mapping[mapping['dist_to_tss'].abs() <= max_dist].copy()

    # Calculate distance weight (exponential decay)
    decay = p2g_cfg.get('decay_distance', 50000)
    mapping['weight'] = np.exp(-mapping['dist_to_tss'].abs() / decay)

    # Boost weight for promoter peaks
    promoter_mask = mapping['annotation'].str.contains('Promoter', case=False, na=False)
    mapping.loc[promoter_mask, 'weight'] = 1.0

    # Clean gene names
    mapping['gene'] = mapping['gene'].astype(str).str.strip()
    mapping = mapping[mapping['gene'] != 'nan'].copy()

    print(f"  Created {len(mapping)} peak-gene links for {mapping['gene'].nunique()} unique genes")

    return mapping


def compute_tf_activity(bindetect: pd.DataFrame, contrast_cfg: dict, scoring_cfg: dict) -> pd.DataFrame:
    """
    Compute TF activity scores from BINDetect results for a given contrast.

    The bindetect file has columns like:
    - "IL10R_IL10"_mean_score, "IL10R_ctrl"_mean_score, etc.
    - "IL10R_IL10"_"IL10R_ctrl"_change, "IL10R_IL10"_"IL10R_ctrl"_pvalue
    - "MPL-JAK2VF_IL10"_"MPL-JAK2VF_ctrl"_change, etc.

    Returns DataFrame with columns: tf, tf_activity, tf_pvalue, direction
    """
    print("\nComputing TF activity scores...")

    all_cols = bindetect.columns.tolist()

    # Parse available groups from mean_score columns
    # Columns look like: "GroupName"_mean_score (with Unicode smart quotes!)
    score_cols = [c for c in all_cols if c.endswith('_mean_score')]
    # Remove Unicode smart quotes (U+201C and U+201D) and regular quotes from group names
    available_groups = []
    for c in score_cols:
        group = c.replace('_mean_score', '')
        group = group.replace('\u201c', '').replace('\u201d', '')  # Unicode smart quotes
        group = group.replace('"', '')  # Regular quotes
        available_groups.append(group)
    print(f"  Available groups in bindetect: {available_groups}")

    # For NvT contrast: Compare MPL-JAK2VF + IL10R-JAK2VF (transformed) vs Parental + IL10R (normal)
    # For CvNC contrast: Compare MPL-JAK2VF vs IL10R-JAK2VF

    contrast_name = contrast_cfg.get('name', '')

    # Define group mappings based on contrast
    if 'Normal' in contrast_name or 'Transformed' in contrast_name or contrast_name == 'Normal_vs_Transformed':
        # NvT: Transformed (JAK2VF samples) vs Normal (Parental, IL10R)
        # Note: bindetect doesn't distinguish CI vs non-CI, so we use all JAK2VF as "transformed"
        case_groups = ['MPL-JAK2VF_IL10', 'MPL-JAK2VF_ctrl', 'IL10R-JAK2VF_IL10', 'IL10R-JAK2VF_ctrl']
        ctrl_groups = ['Parental_IL10', 'Parental_ctrl', 'IL10R_IL10', 'IL10R_ctrl']
        print("  Using NvT contrast: JAK2VF samples vs Parental/IL10R")

    elif 'Canonical' in contrast_name or 'NonCanonical' in contrast_name or contrast_name == 'Canonical_vs_NonCanonical':
        # CvNC: MPL-JAK2VF (non-canonical) vs IL10R-JAK2VF (canonical)
        case_groups = ['MPL-JAK2VF_IL10', 'MPL-JAK2VF_ctrl']  # non-canonical
        ctrl_groups = ['IL10R-JAK2VF_IL10', 'IL10R-JAK2VF_ctrl']  # canonical
        print("  Using CvNC contrast: MPL-JAK2VF vs IL10R-JAK2VF")

    else:
        # Fall back to using all samples
        print(f"  Unknown contrast '{contrast_name}', using all samples")
        case_groups = available_groups[:len(available_groups)//2]
        ctrl_groups = available_groups[len(available_groups)//2:]

    # Filter to available groups
    case_groups = [g for g in case_groups if g in available_groups]
    ctrl_groups = [g for g in ctrl_groups if g in available_groups]

    print(f"  Case groups: {case_groups}")
    print(f"  Control groups: {ctrl_groups}")

    if not case_groups or not ctrl_groups:
        print("  WARNING: No valid groups found, using all TFs with mean activity")
        tf_activity = bindetect[['name']].copy()
        tf_activity['tf_activity'] = bindetect[score_cols].mean(axis=1)
        tf_activity['tf_pvalue'] = 0.01
        tf_activity['direction'] = 1
        tf_activity = tf_activity.rename(columns={'name': 'tf'})
        return tf_activity

    # Look for direct pairwise comparison columns
    # Format: "group1"_"group2"_change
    change_cols = []
    pvalue_cols = []

    for case_g in case_groups:
        for ctrl_g in ctrl_groups:
            # Try both orderings with quotes
            patterns = [
                f'"{case_g}"_"{ctrl_g}"_change',
                f'"{ctrl_g}"_"{case_g}"_change',
                f'{case_g}_{ctrl_g}_change',
                f'{ctrl_g}_{case_g}_change',
            ]
            for p in patterns:
                if p in all_cols:
                    change_cols.append(p)
                    pval_p = p.replace('_change', '_pvalue')
                    if pval_p in all_cols:
                        pvalue_cols.append(pval_p)
                    break

    print(f"  Found {len(change_cols)} pairwise comparison columns")

    if change_cols:
        # Use direct comparisons
        tf_activity = bindetect[['name']].copy()
        tf_activity['tf_activity'] = bindetect[change_cols].mean(axis=1)
        if pvalue_cols:
            tf_activity['tf_pvalue'] = bindetect[pvalue_cols].min(axis=1)
        else:
            tf_activity['tf_pvalue'] = 0.01
    else:
        # Compute from mean scores
        print("  No direct comparisons, computing from mean scores...")

        # Handle quoted column names - the actual column format uses Unicode smart quotes: "GroupName"_mean_score
        case_score_cols = []
        ctrl_score_cols = []

        for g in case_groups:
            # Try Unicode smart quotes first, then regular quotes, then no quotes
            col_patterns = [
                f'\u201c{g}\u201d_mean_score',  # Unicode smart quotes
                f'"{g}"_mean_score',             # Regular quotes
                f'{g}_mean_score'                # No quotes
            ]
            for col_name in col_patterns:
                if col_name in all_cols:
                    case_score_cols.append(col_name)
                    break

        for g in ctrl_groups:
            col_patterns = [
                f'\u201c{g}\u201d_mean_score',
                f'"{g}"_mean_score',
                f'{g}_mean_score'
            ]
            for col_name in col_patterns:
                if col_name in all_cols:
                    ctrl_score_cols.append(col_name)
                    break

        print(f"    Case score columns: {case_score_cols}")
        print(f"    Control score columns: {ctrl_score_cols}")

        if case_score_cols and ctrl_score_cols:
            tf_activity = bindetect[['name']].copy()
            tf_activity['case_mean'] = bindetect[case_score_cols].mean(axis=1)
            tf_activity['ctrl_mean'] = bindetect[ctrl_score_cols].mean(axis=1)
            tf_activity['tf_activity'] = tf_activity['case_mean'] - tf_activity['ctrl_mean']
            tf_activity['tf_pvalue'] = 0.01
        else:
            print("  ERROR: Could not find score columns")
            # Last resort: use all mean scores
            tf_activity = bindetect[['name']].copy()
            tf_activity['tf_activity'] = bindetect[score_cols].std(axis=1)  # use variance as activity proxy
            tf_activity['tf_pvalue'] = 0.01

    tf_activity['direction'] = np.sign(tf_activity['tf_activity'])
    tf_activity = tf_activity.rename(columns={'name': 'tf'})

    # Keep all TFs (don't filter by significance here, filter edges later)
    min_change = scoring_cfg.get('tf_change_min', 0.01)
    significant = tf_activity['tf_activity'].abs() >= min_change
    print(f"  {significant.sum()} TFs have |activity| >= {min_change}")
    print(f"  Top TFs by activity:")
    top = tf_activity.nlargest(5, 'tf_activity')[['tf', 'tf_activity']]
    for _, row in top.iterrows():
        print(f"    {row['tf']}: {row['tf_activity']:.4f}")
    bottom = tf_activity.nsmallest(5, 'tf_activity')[['tf', 'tf_activity']]
    for _, row in bottom.iterrows():
        print(f"    {row['tf']}: {row['tf_activity']:.4f}")

    return tf_activity


def compute_gene_accessibility(atac_counts: pd.DataFrame, metadata: pd.DataFrame,
                               contrast_cfg: dict, scoring_cfg: dict) -> pd.DataFrame:
    """
    Compute gene-level accessibility changes for a given contrast.

    Returns DataFrame with columns: gene, atac_log2fc
    """
    print("\nComputing gene accessibility changes...")

    rna_filter = contrast_cfg.get('rna_filter', {})

    # Map ATAC sample names to conditions
    # ATAC samples are like: PO1, PO2, PIL10, 6O1, etc.
    # Need to map these to metadata conditions

    # For now, use a simple approach based on sample naming
    atac_samples = atac_counts.columns.tolist()

    # Group samples by condition
    # Extract condition info from sample names
    case_samples = []
    ctrl_samples = []

    condition_col = rna_filter.get('column', 'Condition')
    case_val = rna_filter.get('case', 'Transformed')
    ctrl_val = rna_filter.get('control', 'Normal')

    # Map ATAC sample patterns to conditions
    # CI samples = Transformed, non-CI = Normal
    for s in atac_samples:
        if 'CI' in s.upper():
            case_samples.append(s)
        else:
            ctrl_samples.append(s)

    if not case_samples or not ctrl_samples:
        print("  Warning: Could not identify case/control samples from names")
        # Return empty with all genes
        return pd.DataFrame({'gene': atac_counts.index, 'atac_log2fc': 0.0})

    print(f"  Case samples (n={len(case_samples)}): {case_samples[:3]}...")
    print(f"  Control samples (n={len(ctrl_samples)}): {ctrl_samples[:3]}...")

    # Compute log2FC
    case_mean = atac_counts[case_samples].mean(axis=1) + 1  # pseudocount
    ctrl_mean = atac_counts[ctrl_samples].mean(axis=1) + 1
    log2fc = np.log2(case_mean / ctrl_mean)

    result = pd.DataFrame({
        'gene': atac_counts.index,
        'atac_log2fc': log2fc.values
    })

    print(f"  Computed accessibility for {len(result)} genes")

    return result


def compute_rna_de(rna_counts: pd.DataFrame, metadata: pd.DataFrame,
                   contrast_cfg: dict, scoring_cfg: dict) -> pd.DataFrame:
    """
    Compute RNA differential expression for a given contrast.
    Uses a simple fold-change approach (DESeq2 would be better for production).

    Returns DataFrame with columns: gene, rna_log2fc, rna_pvalue
    """
    print("\nComputing RNA differential expression...")

    rna_filter = contrast_cfg.get('rna_filter', {})

    # Get column info
    condition_col = rna_filter.get('column', 'Condition')
    case_val = rna_filter.get('case', 'Transformed')
    ctrl_val = rna_filter.get('control', 'Normal')

    # Handle pre-filter if present
    pre_filter = rna_filter.get('pre_filter', {})
    if pre_filter:
        pre_col = pre_filter.get('column')
        pre_val = pre_filter.get('value')
        if pre_col and pre_val:
            metadata = metadata[metadata[pre_col] == pre_val].copy()
            print(f"  Pre-filtered to {pre_col} == {pre_val}: {len(metadata)} samples")

    # Map sample names between RNA and metadata
    # RNA columns are like: rP_O_1, r6_O_1, r8_CI_O_1, etc.
    # Metadata Probe column is like: P_O_1, 6_O_1, 8_CI_O_1, etc.

    rna_samples = [c for c in rna_counts.columns if c not in ['Description', 'Probe']]

    # Create mapping from RNA column to metadata
    rna_to_meta = {}
    for rna_col in rna_samples:
        # Remove 'r' prefix if present
        meta_name = rna_col[1:] if rna_col.startswith('r') else rna_col
        rna_to_meta[rna_col] = meta_name

    # Find case and control samples
    case_samples = []
    ctrl_samples = []

    for rna_col, meta_name in rna_to_meta.items():
        meta_row = metadata[metadata['Probe'] == meta_name]
        if len(meta_row) == 0:
            continue

        condition = meta_row[condition_col].values[0]
        if condition == case_val:
            case_samples.append(rna_col)
        elif condition == ctrl_val:
            ctrl_samples.append(rna_col)

    print(f"  Case samples ({case_val}, n={len(case_samples)})")
    print(f"  Control samples ({ctrl_val}, n={len(ctrl_samples)})")

    if not case_samples or not ctrl_samples:
        print("  Warning: Could not find samples for contrast")
        return pd.DataFrame({'gene': rna_counts.index, 'rna_log2fc': 0.0, 'rna_pvalue': 1.0})

    # Compute log2FC and simple t-test p-value
    case_data = rna_counts[case_samples]
    ctrl_data = rna_counts[ctrl_samples]

    case_mean = case_data.mean(axis=1) + 1  # pseudocount
    ctrl_mean = ctrl_data.mean(axis=1) + 1
    log2fc = np.log2(case_mean / ctrl_mean)

    # Simple t-test (for a proper analysis, use DESeq2)
    pvalues = []
    for gene in rna_counts.index:
        case_vals = case_data.loc[gene].values
        ctrl_vals = ctrl_data.loc[gene].values

        # Skip if all zeros
        if case_vals.sum() == 0 and ctrl_vals.sum() == 0:
            pvalues.append(1.0)
            continue

        try:
            _, pval = stats.ttest_ind(case_vals, ctrl_vals, equal_var=False)
            pvalues.append(pval if not np.isnan(pval) else 1.0)
        except:
            pvalues.append(1.0)

    result = pd.DataFrame({
        'gene': rna_counts.index,
        'rna_log2fc': log2fc.values,
        'rna_pvalue': pvalues
    })

    # Apply thresholds
    lfc_min = scoring_cfg.get('rna_log2fc_min', 0.5)
    pval_max = scoring_cfg.get('rna_padj_max', 0.05)

    sig_up = (result['rna_log2fc'] >= lfc_min) & (result['rna_pvalue'] <= pval_max)
    sig_down = (result['rna_log2fc'] <= -lfc_min) & (result['rna_pvalue'] <= pval_max)

    print(f"  Significant genes: {sig_up.sum()} up, {sig_down.sum()} down")

    return result


def build_tf_gene_edges(tf_activity: pd.DataFrame, peak2gene: pd.DataFrame,
                        gene_access: pd.DataFrame, rna_de: pd.DataFrame,
                        cfg: dict) -> pd.DataFrame:
    """
    Build TF→gene edges by integrating TF binding, accessibility, and expression.

    EdgeScore = TF_activity * ATAC_accessibility * weight
    Filtered by RNA expression evidence.
    """
    print("\nBuilding TF→gene edges...")

    scoring_cfg = cfg['scoring']

    # Check if we have TFs
    if len(tf_activity) == 0:
        print("  ERROR: No TF activity data!")
        return pd.DataFrame(columns=['tf', 'gene', 'edge_score', 'n_peaks',
                                     'tf_activity', 'atac_log2fc', 'rna_log2fc', 'rna_pvalue'])

    # Get unique genes with their best peak weights
    gene_weights = peak2gene.groupby('gene').agg({
        'weight': 'max',  # best peak weight
        'peak_id': 'count'  # number of peaks
    }).reset_index()
    gene_weights.columns = ['gene', 'peak_weight', 'n_peaks']

    tfs = tf_activity['tf'].unique()
    genes = gene_weights['gene'].unique()

    print(f"  Computing edges for {len(tfs)} TFs x {len(genes)} genes...")

    # Merge gene-level data
    gene_data = gene_weights.merge(gene_access, on='gene', how='left')
    gene_data = gene_data.merge(rna_de, on='gene', how='left')
    gene_data = gene_data.fillna({'atac_log2fc': 0, 'rna_log2fc': 0, 'rna_pvalue': 1})

    # FILTER to differentially expressed genes FIRST (reduces computation)
    rna_lfc_min = scoring_cfg.get('rna_log2fc_min', 0.5)
    de_genes = gene_data[gene_data['rna_log2fc'].abs() >= rna_lfc_min].copy()
    print(f"  DE genes (|log2FC| >= {rna_lfc_min}): {len(de_genes)}")

    if len(de_genes) == 0:
        print("  WARNING: No differentially expressed genes found!")
        # Lower threshold and try again
        rna_lfc_min = 0.25
        de_genes = gene_data[gene_data['rna_log2fc'].abs() >= rna_lfc_min].copy()
        print(f"  Trying lower threshold (|log2FC| >= {rna_lfc_min}): {len(de_genes)} genes")

    if len(de_genes) == 0:
        de_genes = gene_data.head(1000).copy()  # just take top genes
        print(f"  Using top 1000 genes as fallback")

    # FILTER to active TFs (reduces computation)
    tf_min = scoring_cfg.get('tf_change_min', 0.01)
    active_tfs = tf_activity[tf_activity['tf_activity'].abs() >= tf_min].copy()
    print(f"  Active TFs (|activity| >= {tf_min}): {len(active_tfs)}")

    if len(active_tfs) == 0:
        print("  WARNING: No active TFs found, using top 100 by variance")
        active_tfs = tf_activity.nlargest(100, 'tf_activity')

    # Create edges list
    edges = []

    # For each TF, compute edge scores for DE genes
    for _, tf_row in active_tfs.iterrows():
        tf = tf_row['tf']
        tf_act = tf_row['tf_activity']

        # Compute edge scores for DE genes
        gene_edges = de_genes.copy()
        gene_edges['tf'] = tf
        gene_edges['tf_activity'] = tf_act

        # Edge score = TF activity * peak weight * (1 + |accessibility change|)
        # This ensures we get edges even if ATAC change is small
        gene_edges['edge_score'] = (
            gene_edges['tf_activity'] *
            gene_edges['peak_weight'] *
            (1 + gene_edges['atac_log2fc'].abs().clip(0, 2))
        )

        edges.append(gene_edges[['tf', 'gene', 'edge_score', 'n_peaks',
                                  'tf_activity', 'atac_log2fc', 'rna_log2fc', 'rna_pvalue']])

    if not edges:
        print("  ERROR: No edges could be computed!")
        return pd.DataFrame(columns=['tf', 'gene', 'edge_score', 'n_peaks',
                                     'tf_activity', 'atac_log2fc', 'rna_log2fc', 'rna_pvalue'])

    edges_df = pd.concat(edges, ignore_index=True)

    # Filter edges by edge score
    print(f"  Raw edges: {len(edges_df)}")
    edge_min = scoring_cfg.get('edge_score_min', 0.01)
    edges_df = edges_df[edges_df['edge_score'].abs() >= edge_min]
    print(f"  After edge score filter (|score| >= {edge_min}): {len(edges_df)}")

    if len(edges_df) == 0:
        print("  WARNING: No edges passed filter, using top edges")
        edges_df = pd.concat(edges, ignore_index=True)
        edges_df = edges_df.nlargest(10000, 'edge_score')

    # Normalize edge scores
    if cfg['scoring'].get('normalize_edges', True):
        max_score = edges_df['edge_score'].abs().max()
        if max_score > 0:
            edges_df['edge_score_norm'] = edges_df['edge_score'] / max_score
        else:
            edges_df['edge_score_norm'] = 0

    print(f"  Final edges: {len(edges_df)} ({edges_df['tf'].nunique()} TFs, {edges_df['gene'].nunique()} genes)")

    return edges_df


def build_nodes_table(edges: pd.DataFrame, tf_activity: pd.DataFrame,
                      rna_de: pd.DataFrame) -> pd.DataFrame:
    """Build nodes table with attributes for network visualization."""
    print("\nBuilding nodes table...")

    # TF nodes
    tf_nodes = tf_activity[['tf', 'tf_activity', 'direction']].copy()
    tf_nodes['node_type'] = 'TF'
    tf_nodes = tf_nodes.rename(columns={'tf': 'node'})

    # Add degree (number of target genes)
    tf_degree = edges.groupby('tf').size().reset_index(name='degree')
    tf_nodes = tf_nodes.merge(tf_degree.rename(columns={'tf': 'node'}), on='node', how='left')
    tf_nodes['degree'] = tf_nodes['degree'].fillna(0).astype(int)

    # Gene nodes (only genes that appear in edges)
    edge_genes = edges['gene'].unique()
    gene_nodes = rna_de[rna_de['gene'].isin(edge_genes)][['gene', 'rna_log2fc']].copy()
    gene_nodes['node_type'] = 'gene'
    gene_nodes = gene_nodes.rename(columns={'gene': 'node', 'rna_log2fc': 'expression'})
    gene_nodes['tf_activity'] = np.nan
    gene_nodes['direction'] = np.sign(gene_nodes['expression'])

    # Add degree (number of regulating TFs)
    gene_degree = edges.groupby('gene').size().reset_index(name='degree')
    gene_nodes = gene_nodes.merge(gene_degree.rename(columns={'gene': 'node'}), on='node', how='left')

    # Combine
    nodes = pd.concat([tf_nodes, gene_nodes], ignore_index=True)

    print(f"  Built {len(nodes)} nodes ({(nodes['node_type']=='TF').sum()} TFs, {(nodes['node_type']=='gene').sum()} genes)")

    return nodes


def main():
    parser = argparse.ArgumentParser(description='Build TF-Gene Regulatory Network')
    parser.add_argument('--config', type=str, default='./config/config.yaml', help='Config file path')
    parser.add_argument('--contrast', type=str, required=True, help='Contrast name (e.g., NvT, CvNC)')
    parser.add_argument('--output', type=str, help='Output directory (overrides config)')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    # Load config
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Error: Config file not found: {config_path}")
        print("Please create a config.yaml file. See config/config_template.yaml for an example.")
        sys.exit(1)

    cfg = load_config(config_path)
    base_dir = config_path.parent

    # Check contrast exists
    if args.contrast not in cfg['contrasts']:
        print(f"Error: Contrast '{args.contrast}' not found in config")
        print(f"Available contrasts: {list(cfg['contrasts'].keys())}")
        sys.exit(1)

    contrast_cfg = cfg['contrasts'][args.contrast]
    print(f"\n{'='*60}")
    print(f"Building TF-Gene Network for: {contrast_cfg['name']}")
    print(f"{'='*60}")
    print(f"Description: {contrast_cfg.get('description', '')}")

    # Setup output directory
    out_dir = Path(args.output) if args.output else base_dir / cfg['output']['base_dir'] / args.contrast
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput directory: {out_dir}")

    # Load data
    print("\n" + "-"*40)
    print("LOADING DATA")
    print("-"*40)

    bindetect = load_bindetect_results(resolve_path(cfg['data']['bindetect_results'], base_dir))
    peak_annot = load_peak_annotations(resolve_path(cfg['data']['atac_peaks_annotated'], base_dir))
    atac_counts = load_atac_gene_counts(resolve_path(cfg['data']['atac_gene_counts'], base_dir))
    rna_counts = load_rna_counts(resolve_path(cfg['data']['rna_counts'], base_dir))
    metadata = load_metadata(resolve_path(cfg['data']['metadata'], base_dir))

    # Create peak-to-gene mapping
    print("\n" + "-"*40)
    print("PROCESSING")
    print("-"*40)

    peak2gene = create_peak2gene_mapping(peak_annot, cfg)

    # Compute TF activity
    tf_activity = compute_tf_activity(bindetect, contrast_cfg, cfg['scoring'])

    # Compute gene accessibility
    gene_access = compute_gene_accessibility(atac_counts, metadata, contrast_cfg, cfg['scoring'])

    # Compute RNA DE
    rna_de = compute_rna_de(rna_counts, metadata, contrast_cfg, cfg['scoring'])

    # Build edges
    edges = build_tf_gene_edges(tf_activity, peak2gene, gene_access, rna_de, cfg)

    # Build nodes
    nodes = build_nodes_table(edges, tf_activity, rna_de)

    # Save outputs
    print("\n" + "-"*40)
    print("SAVING OUTPUTS")
    print("-"*40)

    tables_dir = out_dir / cfg['output']['tables_dir']
    tables_dir.mkdir(exist_ok=True)

    edges_file = tables_dir / f"tf_gene_edges_{args.contrast}.tsv"
    nodes_file = tables_dir / f"tf_gene_nodes_{args.contrast}.tsv"
    peak2gene_file = tables_dir / f"peak2gene_{args.contrast}.tsv"
    tf_activity_file = tables_dir / f"tf_activity_{args.contrast}.tsv"
    rna_de_file = tables_dir / f"rna_de_{args.contrast}.tsv"

    edges.to_csv(edges_file, sep='\t', index=False)
    nodes.to_csv(nodes_file, sep='\t', index=False)
    peak2gene.to_csv(peak2gene_file, sep='\t', index=False)
    tf_activity.to_csv(tf_activity_file, sep='\t', index=False)
    rna_de.to_csv(rna_de_file, sep='\t', index=False)

    print(f"  Edges: {edges_file}")
    print(f"  Nodes: {nodes_file}")
    print(f"  Peak2Gene: {peak2gene_file}")
    print(f"  TF Activity: {tf_activity_file}")
    print(f"  RNA DE: {rna_de_file}")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"  Total edges: {len(edges)}")
    print(f"  Total TFs: {edges['tf'].nunique()}")
    print(f"  Total genes: {edges['gene'].nunique()}")
    print(f"\n  Top 10 TFs by activity:")
    top_tfs = tf_activity.nlargest(10, 'tf_activity')[['tf', 'tf_activity']]
    for _, row in top_tfs.iterrows():
        print(f"    {row['tf']}: {row['tf_activity']:.4f}")

    print(f"\n  Top 10 TFs by number of targets:")
    top_by_targets = edges.groupby('tf').size().nlargest(10)
    for tf, count in top_by_targets.items():
        print(f"    {tf}: {count} targets")

    print("\nDone!")


if __name__ == '__main__':
    main()
