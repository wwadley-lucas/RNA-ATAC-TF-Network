#!/usr/bin/env python3
"""
TF-Focused Regulatory Network Builder
======================================
Builds cascading regulatory networks centered on specific TFs.
Shows primary TF → direct targets, and if targets are TFs themselves,
expands to show their downstream targets (2 layers deep).

Part of RNA-ATAC-TF-Network pipeline.
https://github.com/your-username/RNA-ATAC-TF-Network

Usage:
    python build_tf_focus_network.py --config config.yaml --contrast NvT
    python build_tf_focus_network.py --config config.yaml --tf STAT1

Output:
    - tf_focus_edges_{tf}_{contrast}.tsv: Edge table with cascade info
    - tf_focus_nodes_{tf}_{contrast}.tsv: Node attributes for visualization
"""

__version__ = "1.0.0"

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import yaml


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


def load_pathway_genes(pathway_name: str, species: str = "Mus musculus", cache_dir: Optional[Path] = None) -> Set[str]:
    """
    Load genes from MSigDB pathway using msigdbr via R subprocess.

    Args:
        pathway_name: Name of the pathway (e.g., "HALLMARK_INTERFERON_ALPHA_RESPONSE")
        species: Species name ("Mus musculus" or "Homo sapiens")
        cache_dir: Optional directory to cache pathway gene lists

    Returns:
        Set of gene symbols in the pathway
    """
    # Check for cached file first
    if cache_dir:
        cache_file = cache_dir / f"{pathway_name}_{species.replace(' ', '_')}_genes.txt"
        if cache_file.exists():
            print(f"    Loading cached pathway genes from {cache_file}")
            with open(cache_file, 'r') as f:
                genes = set(line.strip() for line in f if line.strip())
            print(f"    Loaded {len(genes)} genes from {pathway_name}")
            return genes

    # Use R to fetch pathway genes via msigdbr
    r_code = f'''
    suppressPackageStartupMessages({{
        library(msigdbr)
        library(dplyr)
    }})
    genes <- msigdbr(species = "{species}", category = "H") %>%
        filter(gs_name == "{pathway_name}") %>%
        pull(gene_symbol) %>%
        unique()
    cat(genes, sep = "\\n")
    '''

    try:
        result = subprocess.run(
            ['Rscript', '-e', r_code],
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode == 0:
            genes = set(line.strip() for line in result.stdout.split('\n') if line.strip())
            print(f"    Loaded {len(genes)} genes from {pathway_name}")

            # Cache the result
            if cache_dir:
                cache_dir.mkdir(parents=True, exist_ok=True)
                with open(cache_file, 'w') as f:
                    f.write('\n'.join(sorted(genes)))
                print(f"    Cached pathway genes to {cache_file}")

            return genes
        else:
            print(f"    WARNING: Failed to load pathway genes: {result.stderr}")
            return set()

    except subprocess.TimeoutExpired:
        print(f"    WARNING: Timeout loading pathway genes")
        return set()
    except FileNotFoundError:
        print(f"    WARNING: Rscript not found, cannot load pathway genes")
        return set()


def get_tf_list(bindetect_dir: Path) -> Dict[str, str]:
    """
    Scan BINDetect folder to get all TF names and their motif IDs.

    Returns:
        Dict mapping TF name (uppercase) to motif ID
        e.g., {"STAT1": "MA0137.4", "IRF1": "MA0050.4", ...}
    """
    print(f"  Scanning TF folders in {bindetect_dir}...")

    tf_list = {}

    if not bindetect_dir.exists():
        print(f"    WARNING: BINDetect directory not found: {bindetect_dir}")
        return tf_list

    for folder in bindetect_dir.iterdir():
        if folder.is_dir() and '_' in folder.name:
            # Parse folder name: TFname_MotifID (e.g., STAT1_MA0137.4)
            parts = folder.name.rsplit('_', 1)
            if len(parts) == 2:
                tf_name = parts[0].upper()
                motif_id = parts[1]
                # Handle cases like IRF1_MA0050.4 - store both original and uppercase
                tf_list[tf_name] = motif_id
                # Also store with original case
                tf_list[parts[0]] = motif_id

    print(f"    Found {len(tf_list) // 2} unique TFs")
    return tf_list


def get_condition_column_indices() -> Dict[str, int]:
    """
    Map condition names to BED file column indices.

    BED columns 10-17 (0-indexed: 9-16) contain accessibility values for:
    0: Parental_ctrl
    1: Parental_IL10
    2: IL10R_ctrl
    3: IL10R_IL10
    4: IL10R-JAK2VF_ctrl
    5: IL10R-JAK2VF_IL10
    6: MPL-JAK2VF_ctrl
    7: MPL-JAK2VF_IL10
    """
    return {
        "Parental_ctrl": 0,
        "Parental_IL10": 1,
        "IL10R_ctrl": 2,
        "IL10R_IL10": 3,
        "IL10R-JAK2VF_ctrl": 4,
        "IL10R-JAK2VF_IL10": 5,
        "MPL-JAK2VF_ctrl": 6,
        "MPL-JAK2VF_IL10": 7,
    }


def load_tf_bed_file(
    tf_name: str,
    motif_id: str,
    bindetect_dir: Path,
    use_bound_only: bool = True,
    case_conditions: List[str] = None,
    min_binding_score: float = 5.0
) -> pd.DataFrame:
    """
    Load BED file for specific TF.

    Args:
        tf_name: TF name (e.g., "STAT1")
        motif_id: JASPAR motif ID (e.g., "MA0137.4")
        bindetect_dir: Path to BINDetect output directory
        use_bound_only: If True, only load bound sites for case conditions
        case_conditions: List of condition names to consider as "bound"
        min_binding_score: Minimum motif binding score (column 5)

    Returns:
        DataFrame with columns:
            site_id, chr, start, end, binding_score, strand,
            peak_chr, peak_start, peak_end, accessibility_0, ..., accessibility_7
    """
    # Build folder path - try different case variations
    folder_name = f"{tf_name}_{motif_id}"
    folder_path = bindetect_dir / folder_name / "beds"

    if not folder_path.exists():
        # Try original case from folder listing
        for item in bindetect_dir.iterdir():
            if item.is_dir() and motif_id in item.name:
                folder_path = item / "beds"
                folder_name = item.name
                break

    if not folder_path.exists():
        print(f"    WARNING: BED folder not found for {tf_name}_{motif_id}")
        return pd.DataFrame()

    # Determine which file to load
    if use_bound_only and case_conditions:
        # Load bound sites for each case condition and merge
        all_sites = []
        for cond in case_conditions:
            bound_file = folder_path / f'{folder_name}_"{cond}"_bound.bed'
            if bound_file.exists():
                df = _parse_bed_file(bound_file)
                if len(df) > 0:
                    all_sites.append(df)

        if all_sites:
            sites = pd.concat(all_sites, ignore_index=True)
            # Remove duplicates (same site may be bound in multiple conditions)
            sites = sites.drop_duplicates(subset=['chr', 'start', 'end'])
        else:
            # Fall back to all sites
            all_file = folder_path / f"{folder_name}_all.bed"
            sites = _parse_bed_file(all_file) if all_file.exists() else pd.DataFrame()
    else:
        # Load all sites
        all_file = folder_path / f"{folder_name}_all.bed"
        sites = _parse_bed_file(all_file) if all_file.exists() else pd.DataFrame()

    if len(sites) == 0:
        print(f"    WARNING: No binding sites found for {tf_name}")
        return sites

    # Filter by binding score
    if 'binding_score' in sites.columns:
        sites = sites[sites['binding_score'] >= min_binding_score].copy()

    print(f"    Loaded {len(sites)} binding sites for {tf_name}")
    return sites


def _parse_bed_file(bed_path: Path) -> pd.DataFrame:
    """Parse a TOBIAS BED file into a DataFrame."""
    if not bed_path.exists():
        return pd.DataFrame()

    # BED columns:
    # 0: chr, 1: start, 2: end (motif site)
    # 3: TF name
    # 4: binding score (motif affinity)
    # 5: strand
    # 6: peak_chr, 7: peak_start, 8: peak_end
    # 9-16: accessibility values for 8 conditions

    col_names = [
        'chr', 'start', 'end', 'tf_name', 'binding_score', 'strand',
        'peak_chr', 'peak_start', 'peak_end',
        'acc_0', 'acc_1', 'acc_2', 'acc_3', 'acc_4', 'acc_5', 'acc_6', 'acc_7'
    ]

    try:
        df = pd.read_csv(bed_path, sep='\t', header=None, names=col_names)
        df['site_id'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
        return df
    except Exception as e:
        print(f"    Error reading {bed_path}: {e}")
        return pd.DataFrame()


def map_sites_to_genes(
    binding_sites: pd.DataFrame,
    peak2gene_mapping: pd.DataFrame,
    aggregation: str = 'max'
) -> pd.DataFrame:
    """
    Link binding sites to genes via peak coordinates.

    Args:
        binding_sites: DataFrame from load_tf_bed_file
        peak2gene_mapping: DataFrame with peak_id, gene, weight columns
        aggregation: How to aggregate multiple sites per gene ('max' or 'sum')

    Returns:
        DataFrame with columns: gene, binding_score, accessibility scores, n_sites
    """
    if len(binding_sites) == 0:
        return pd.DataFrame()

    # Create peak_id from binding site peak coordinates
    # Note: BED files are 0-based, but peak2gene may be 1-based
    # Try both coordinate systems
    binding_sites = binding_sites.copy()

    # Try 0-based first (BED native)
    binding_sites['peak_id_0based'] = (
        binding_sites['peak_chr'] + ':' +
        binding_sites['peak_start'].astype(str) + '-' +
        binding_sites['peak_end'].astype(str)
    )

    # Try 1-based (start+1, common in annotations)
    binding_sites['peak_id_1based'] = (
        binding_sites['peak_chr'] + ':' +
        (binding_sites['peak_start'] + 1).astype(str) + '-' +
        binding_sites['peak_end'].astype(str)
    )

    # Check which coordinate system matches better
    matches_0based = binding_sites['peak_id_0based'].isin(peak2gene_mapping['peak_id']).sum()
    matches_1based = binding_sites['peak_id_1based'].isin(peak2gene_mapping['peak_id']).sum()

    if matches_1based > matches_0based:
        binding_sites['peak_id'] = binding_sites['peak_id_1based']
        print(f"    Using 1-based coordinates ({matches_1based} matches vs {matches_0based})")
    else:
        binding_sites['peak_id'] = binding_sites['peak_id_0based']
        print(f"    Using 0-based coordinates ({matches_0based} matches vs {matches_1based})")

    # Merge with peak2gene mapping
    merged = binding_sites.merge(
        peak2gene_mapping[['peak_id', 'gene', 'weight']],
        on='peak_id',
        how='inner'
    )

    if len(merged) == 0:
        print("    WARNING: No binding sites matched peak2gene mapping")
        return pd.DataFrame()

    # Aggregate by gene
    acc_cols = [c for c in merged.columns if c.startswith('acc_')]

    if aggregation == 'max':
        agg_dict = {'binding_score': 'max', 'weight': 'max', 'site_id': 'count'}
        for col in acc_cols:
            agg_dict[col] = 'max'
    else:
        agg_dict = {'binding_score': 'sum', 'weight': 'sum', 'site_id': 'count'}
        for col in acc_cols:
            agg_dict[col] = 'mean'

    gene_scores = merged.groupby('gene').agg(agg_dict).reset_index()
    gene_scores = gene_scores.rename(columns={'site_id': 'n_sites'})

    print(f"    Mapped to {len(gene_scores)} unique genes")
    return gene_scores


def compute_contrast_scores(
    tf_gene_scores: pd.DataFrame,
    case_conditions: List[str],
    control_conditions: List[str],
    mode: str = 'differential'
) -> pd.DataFrame:
    """
    Calculate differential or absolute accessibility scores.

    Args:
        tf_gene_scores: DataFrame from map_sites_to_genes
        case_conditions: List of case condition names
        control_conditions: List of control condition names
        mode: 'differential' (case - control) or 'absolute' (case only)

    Returns:
        DataFrame with added accessibility_score column
    """
    if len(tf_gene_scores) == 0:
        return tf_gene_scores

    condition_indices = get_condition_column_indices()

    # Get column names for case and control conditions
    case_cols = [f'acc_{condition_indices[c]}' for c in case_conditions if c in condition_indices]
    ctrl_cols = [f'acc_{condition_indices[c]}' for c in control_conditions if c in condition_indices]

    tf_gene_scores = tf_gene_scores.copy()

    if mode == 'differential':
        case_mean = tf_gene_scores[case_cols].mean(axis=1) if case_cols else 0
        ctrl_mean = tf_gene_scores[ctrl_cols].mean(axis=1) if ctrl_cols else 0
        tf_gene_scores['accessibility_score'] = case_mean - ctrl_mean
    else:
        # Absolute mode - just use case values
        tf_gene_scores['accessibility_score'] = tf_gene_scores[case_cols].mean(axis=1) if case_cols else 0

    return tf_gene_scores


def filter_to_de_genes(
    tf_gene_scores: pd.DataFrame,
    rna_de: pd.DataFrame,
    lfc_threshold: float = 0.5,
    pval_threshold: float = 0.05
) -> pd.DataFrame:
    """
    Filter to differentially expressed genes and merge expression data.

    Args:
        tf_gene_scores: DataFrame from compute_contrast_scores
        rna_de: DataFrame with gene, rna_log2fc, rna_pvalue columns
        lfc_threshold: |log2FC| threshold for DE
        pval_threshold: p-value threshold for DE

    Returns:
        DataFrame with added rna_log2fc column, filtered to DE genes
    """
    if len(tf_gene_scores) == 0:
        return tf_gene_scores

    # Merge with RNA DE data
    merged = tf_gene_scores.merge(
        rna_de[['gene', 'rna_log2fc', 'rna_pvalue']],
        on='gene',
        how='left'
    )

    # Fill NAs
    merged['rna_log2fc'] = merged['rna_log2fc'].fillna(0)
    merged['rna_pvalue'] = merged['rna_pvalue'].fillna(1)

    # Filter to DE genes
    de_mask = (merged['rna_log2fc'].abs() >= lfc_threshold) & (merged['rna_pvalue'] <= pval_threshold)
    filtered = merged[de_mask].copy()

    print(f"    Filtered to {len(filtered)} DE genes (|log2FC| >= {lfc_threshold}, p <= {pval_threshold})")
    return filtered


def get_tf_expression(
    tf_name: str,
    rna_de: pd.DataFrame
) -> Tuple[float, float]:
    """
    Get TF's own expression log2FC and p-value.

    Returns:
        Tuple of (log2FC, pvalue)
    """
    tf_row = rna_de[rna_de['gene'].str.upper() == tf_name.upper()]

    if len(tf_row) > 0:
        return tf_row['rna_log2fc'].values[0], tf_row['rna_pvalue'].values[0]
    return 0.0, 1.0


def build_cascade_network(
    primary_tf: str,
    primary_motif: str,
    tf_list: Dict[str, str],
    peak2gene: pd.DataFrame,
    rna_de: pd.DataFrame,
    bindetect_dir: Path,
    config: dict
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build cascading regulatory network from primary TF.

    Steps:
    1. Get primary TF's DE targets
    2. Identify which targets are TFs (using tf_list lookup)
    3. For each secondary TF, get their DE targets
    4. Build combined edge table with layer info

    Args:
        primary_tf: Name of primary TF (e.g., "STAT1")
        primary_motif: Motif ID (e.g., "MA0137.4")
        tf_list: Dict mapping TF names to motif IDs
        peak2gene: Peak-to-gene mapping DataFrame
        rna_de: RNA DE results DataFrame
        bindetect_dir: Path to BINDetect directory
        config: tf_focus config section

    Returns:
        Tuple of (edges_df, nodes_df)
    """
    print(f"\n  Building cascade network for {primary_tf}...")

    # Config settings
    tf_cfg = config
    cascade_cfg = tf_cfg.get('cascade', {})
    contrast_cfg = tf_cfg.get('contrast', {})

    case_conditions = contrast_cfg.get('case_conditions', [])
    control_conditions = contrast_cfg.get('control_conditions', [])
    accessibility_mode = tf_cfg.get('accessibility_mode', 'differential')
    use_bound_only = tf_cfg.get('use_bound_sites_only', True)
    min_binding_score = tf_cfg.get('min_binding_score', 5.0)
    only_de = tf_cfg.get('only_de_genes', True)
    de_lfc = tf_cfg.get('de_log2fc_threshold', 0.5)
    de_pval = tf_cfg.get('de_pvalue_threshold', 0.05)
    max_targets = tf_cfg.get('max_target_genes', 50)

    cascade_enabled = cascade_cfg.get('enabled', True)
    max_secondary_tfs = cascade_cfg.get('max_secondary_tfs', 5)
    max_targets_per_secondary = cascade_cfg.get('max_targets_per_secondary', 20)
    invert_expression = contrast_cfg.get('invert_expression', False)

    # Pathway filtering settings
    pathway_cfg = tf_cfg.get('pathway_filter', {})
    pathway_enabled = pathway_cfg.get('enabled', False)
    pathway_name = pathway_cfg.get('pathway', None)
    pathway_genes: Set[str] = set()

    if pathway_enabled and pathway_name:
        print(f"  Loading pathway genes: {pathway_name}")
        species = tf_cfg.get('species', 'Mus musculus')
        cache_dir = bindetect_dir.parent / "pathway_cache"
        pathway_genes = load_pathway_genes(pathway_name, species, cache_dir)

    all_edges = []
    all_nodes = {}

    # Invert RNA expression if configured (to match contrast direction)
    if invert_expression:
        print("  Inverting RNA expression values to match contrast direction")
        rna_de = rna_de.copy()
        rna_de['rna_log2fc'] = -rna_de['rna_log2fc']

    # ========== Layer 0: Primary TF ==========
    print(f"  Layer 0: Primary TF {primary_tf}")

    # Get primary TF expression
    primary_lfc, primary_pval = get_tf_expression(primary_tf, rna_de)

    all_nodes[primary_tf] = {
        'node': primary_tf,
        'node_type': 'TF',
        'layer': 0,
        'expression': primary_lfc,
        'expression_pval': primary_pval,
        'is_primary': True
    }

    # ========== Layer 1: Primary TF targets ==========
    print(f"  Layer 1: Loading {primary_tf} targets...")

    # Load primary TF binding sites
    primary_sites = load_tf_bed_file(
        primary_tf, primary_motif, bindetect_dir,
        use_bound_only=use_bound_only,
        case_conditions=case_conditions,
        min_binding_score=min_binding_score
    )

    if len(primary_sites) == 0:
        print(f"    ERROR: No binding sites found for {primary_tf}")
        return pd.DataFrame(), pd.DataFrame()

    # Map to genes
    primary_gene_scores = map_sites_to_genes(primary_sites, peak2gene)

    if len(primary_gene_scores) == 0:
        print(f"    ERROR: No genes mapped for {primary_tf}")
        return pd.DataFrame(), pd.DataFrame()

    # Compute accessibility scores
    primary_gene_scores = compute_contrast_scores(
        primary_gene_scores, case_conditions, control_conditions, accessibility_mode
    )

    # Filter to DE genes if requested
    if only_de:
        primary_targets = filter_to_de_genes(primary_gene_scores, rna_de, de_lfc, de_pval)
    else:
        primary_targets = primary_gene_scores.merge(
            rna_de[['gene', 'rna_log2fc', 'rna_pvalue']], on='gene', how='left'
        )
        primary_targets['rna_log2fc'] = primary_targets['rna_log2fc'].fillna(0)

    # Filter to pathway genes if enabled
    if pathway_enabled and pathway_genes:
        pre_filter_count = len(primary_targets)
        # Keep genes that are in the pathway (case-insensitive match)
        pathway_genes_upper = {g.upper() for g in pathway_genes}
        primary_targets = primary_targets[
            primary_targets['gene'].str.upper().isin(pathway_genes_upper)
        ].copy()
        print(f"    Filtered to pathway genes: {pre_filter_count} -> {len(primary_targets)}")

    if len(primary_targets) == 0:
        print(f"    WARNING: No DE targets found for {primary_tf}, using top genes")
        primary_targets = primary_gene_scores.nlargest(max_targets, 'binding_score')
        primary_targets = primary_targets.merge(
            rna_de[['gene', 'rna_log2fc', 'rna_pvalue']], on='gene', how='left'
        )
        primary_targets['rna_log2fc'] = primary_targets['rna_log2fc'].fillna(0)

    # Rank by combined score
    primary_targets['combined_score'] = (
        primary_targets['binding_score'] *
        (1 + primary_targets['accessibility_score'].abs()) *
        (1 + primary_targets['rna_log2fc'].abs())
    )

    # Identify which targets are TFs (before limiting)
    primary_targets['is_tf'] = primary_targets['gene'].apply(
        lambda g: g.upper() in tf_list and g.upper() != primary_tf.upper()
    )

    # Separate TFs and non-TFs to ensure TFs get included
    tf_targets = primary_targets[primary_targets['is_tf']].copy()
    non_tf_targets = primary_targets[~primary_targets['is_tf']].copy()

    # Take top TFs (by combined score) and top non-TFs
    max_secondary = cascade_cfg.get('max_secondary_tfs', 5) if cascade_enabled else 0
    n_tf_slots = min(max_secondary, len(tf_targets))
    n_gene_slots = max_targets - n_tf_slots

    top_tfs = tf_targets.nlargest(n_tf_slots, 'combined_score') if n_tf_slots > 0 else pd.DataFrame()
    top_genes = non_tf_targets.nlargest(n_gene_slots, 'combined_score')

    # Combine and sort
    primary_targets = pd.concat([top_tfs, top_genes], ignore_index=True)
    primary_targets = primary_targets.sort_values('combined_score', ascending=False)

    # Identify secondary TFs for cascade
    secondary_tfs = []
    for _, row in top_tfs.iterrows():
        gene = row['gene']
        gene_upper = gene.upper()
        if gene_upper in tf_list:
            secondary_tfs.append({
                'name': gene,
                'motif_id': tf_list[gene_upper]
            })

    print(f"    Found {len(secondary_tfs)} secondary TFs among {len(primary_targets)} targets")

    # Create Layer 1 edges and nodes
    for _, row in primary_targets.iterrows():
        gene = row['gene']
        gene_upper = gene.upper()

        # Determine if this gene is a TF
        is_tf = gene_upper in tf_list and gene_upper != primary_tf.upper()

        # Edge
        all_edges.append({
            'source': primary_tf,
            'target': gene,
            'binding_score': row['binding_score'],
            'accessibility_diff': row['accessibility_score'],
            'n_sites': row.get('n_sites', 1),
            'layer': 1,
            'source_type': 'TF',
            'target_type': 'TF' if is_tf else 'gene'
        })

        # Node
        if gene not in all_nodes:
            all_nodes[gene] = {
                'node': gene,
                'node_type': 'TF' if is_tf else 'gene',
                'layer': 1,
                'expression': row['rna_log2fc'],
                'expression_pval': row.get('rna_pvalue', 1.0),
                'is_primary': False
            }

    # ========== Layer 2: Secondary TF targets ==========
    if cascade_enabled and secondary_tfs:
        print(f"\n  Layer 2: Expanding secondary TFs...")

        # Limit number of secondary TFs to expand
        secondary_tfs = secondary_tfs[:max_secondary_tfs]

        for sec_tf in secondary_tfs:
            sec_name = sec_tf['name']
            sec_motif = sec_tf['motif_id']

            print(f"    Expanding {sec_name}...")

            # Load secondary TF binding sites
            sec_sites = load_tf_bed_file(
                sec_name, sec_motif, bindetect_dir,
                use_bound_only=use_bound_only,
                case_conditions=case_conditions,
                min_binding_score=min_binding_score
            )

            if len(sec_sites) == 0:
                continue

            # Map to genes
            sec_gene_scores = map_sites_to_genes(sec_sites, peak2gene)

            if len(sec_gene_scores) == 0:
                continue

            # Compute accessibility scores
            sec_gene_scores = compute_contrast_scores(
                sec_gene_scores, case_conditions, control_conditions, accessibility_mode
            )

            # Filter to DE genes
            if only_de:
                sec_targets = filter_to_de_genes(sec_gene_scores, rna_de, de_lfc, de_pval)
            else:
                sec_targets = sec_gene_scores.merge(
                    rna_de[['gene', 'rna_log2fc', 'rna_pvalue']], on='gene', how='left'
                )
                sec_targets['rna_log2fc'] = sec_targets['rna_log2fc'].fillna(0)

            # Filter to pathway genes if enabled
            if pathway_enabled and pathway_genes:
                pathway_genes_upper = {g.upper() for g in pathway_genes}
                sec_targets = sec_targets[
                    sec_targets['gene'].str.upper().isin(pathway_genes_upper)
                ].copy()

            # Limit targets per secondary TF
            if len(sec_targets) > max_targets_per_secondary:
                sec_targets['combined_score'] = (
                    sec_targets['binding_score'] *
                    (1 + sec_targets['accessibility_score'].abs()) *
                    (1 + sec_targets['rna_log2fc'].abs())
                )
                sec_targets = sec_targets.nlargest(max_targets_per_secondary, 'combined_score')

            # Create Layer 2 edges and nodes
            for _, row in sec_targets.iterrows():
                target_gene = row['gene']

                # Skip if target is the primary TF or already a layer 1 node
                if target_gene.upper() == primary_tf.upper():
                    continue

                # Edge
                all_edges.append({
                    'source': sec_name,
                    'target': target_gene,
                    'binding_score': row['binding_score'],
                    'accessibility_diff': row['accessibility_score'],
                    'n_sites': row.get('n_sites', 1),
                    'layer': 2,
                    'source_type': 'TF',
                    'target_type': 'gene'
                })

                # Node (only add if not already present)
                if target_gene not in all_nodes:
                    all_nodes[target_gene] = {
                        'node': target_gene,
                        'node_type': 'gene',
                        'layer': 2,
                        'expression': row['rna_log2fc'],
                        'expression_pval': row.get('rna_pvalue', 1.0),
                        'is_primary': False
                    }

    # Convert to DataFrames
    edges_df = pd.DataFrame(all_edges) if all_edges else pd.DataFrame()
    nodes_df = pd.DataFrame(list(all_nodes.values())) if all_nodes else pd.DataFrame()

    # Add scaled scores for visualization
    if len(edges_df) > 0:
        # Scale binding score to 0-1
        max_binding = edges_df['binding_score'].max()
        if max_binding > 0:
            edges_df['binding_score_scaled'] = edges_df['binding_score'] / max_binding
        else:
            edges_df['binding_score_scaled'] = 0.5

        # Scale accessibility diff (handle negative values)
        max_acc = edges_df['accessibility_diff'].abs().max()
        if max_acc > 0:
            edges_df['accessibility_diff_scaled'] = edges_df['accessibility_diff'] / max_acc
        else:
            edges_df['accessibility_diff_scaled'] = 0

    # Add degree info to nodes
    if len(edges_df) > 0 and len(nodes_df) > 0:
        out_degree = edges_df.groupby('source').size().reset_index(name='out_degree')
        in_degree = edges_df.groupby('target').size().reset_index(name='in_degree')

        nodes_df = nodes_df.merge(out_degree.rename(columns={'source': 'node'}), on='node', how='left')
        nodes_df = nodes_df.merge(in_degree.rename(columns={'target': 'node'}), on='node', how='left')
        nodes_df['out_degree'] = nodes_df['out_degree'].fillna(0).astype(int)
        nodes_df['in_degree'] = nodes_df['in_degree'].fillna(0).astype(int)
        nodes_df['degree'] = nodes_df['out_degree'] + nodes_df['in_degree']

    # Summary
    print(f"\n  Network summary:")
    print(f"    Total nodes: {len(nodes_df)}")
    print(f"    Total edges: {len(edges_df)}")

    if len(nodes_df) > 0:
        print(f"    Layer 0 (primary TF): {len(nodes_df[nodes_df['layer'] == 0])}")
        print(f"    Layer 1 (direct targets): {len(nodes_df[nodes_df['layer'] == 1])}")
        print(f"    Layer 2 (secondary targets): {len(nodes_df[nodes_df['layer'] == 2])}")
        print(f"    TFs: {len(nodes_df[nodes_df['node_type'] == 'TF'])}")
        print(f"    Genes: {len(nodes_df[nodes_df['node_type'] == 'gene'])}")

    return edges_df, nodes_df


def main():
    parser = argparse.ArgumentParser(description='Build TF-Focused Regulatory Network')
    parser.add_argument('--config', type=str, default='./config/config.yaml', help='Config file path')
    parser.add_argument('--contrast', type=str, default=None, help='Contrast name (uses NvT for RNA DE if not specified)')
    parser.add_argument('--tf', type=str, default=None, help='Specific TF to focus on (overrides config)')
    parser.add_argument('--output', type=str, default=None, help='Output directory (overrides config)')
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

    # Check if tf_focus is enabled
    tf_cfg = cfg.get('tf_focus', {})
    if not tf_cfg.get('enabled', False) and not args.tf:
        print("TF-focus mode is disabled in config and no --tf specified. Set tf_focus.enabled=true or use --tf.")
        sys.exit(0)

    print("\n" + "="*60)
    print("TF-Focused Regulatory Network Builder")
    print("="*60)

    # Get BINDetect directory
    bindetect_dir = Path(tf_cfg.get('bindetect_bed_dir', ''))
    if not bindetect_dir.exists():
        print(f"ERROR: BINDetect directory not found: {bindetect_dir}")
        sys.exit(1)

    print(f"\nBINDetect directory: {bindetect_dir}")

    # Get list of all TFs
    tf_list = get_tf_list(bindetect_dir)

    # Get TFs to process
    if args.tf:
        # Override from command line
        tf_upper = args.tf.upper()
        if tf_upper in tf_list:
            tfs_to_process = [{'name': args.tf, 'motif_id': tf_list[tf_upper]}]
        else:
            print(f"ERROR: TF '{args.tf}' not found in BINDetect directory")
            print(f"Available TFs (first 20): {list(tf_list.keys())[:20]}")
            sys.exit(1)
    else:
        tfs_to_process = tf_cfg.get('tfs', [])

    if not tfs_to_process:
        print("ERROR: No TFs specified in config or command line")
        sys.exit(1)

    # Load RNA DE data (use existing pipeline output)
    contrast = args.contrast or 'NvT'
    rna_de_path = base_dir / cfg['output']['base_dir'] / contrast / cfg['output']['tables_dir'] / f"rna_de_{contrast}.tsv"

    if not rna_de_path.exists():
        print(f"WARNING: RNA DE file not found: {rna_de_path}")
        print("  Attempting to compute RNA DE...")
        # Fall back to loading from raw data
        rna_counts_path = resolve_path(cfg['data']['rna_counts'], base_dir)
        metadata_path = resolve_path(cfg['data']['metadata'], base_dir)

        if rna_counts_path.exists() and metadata_path.exists():
            # Import from main pipeline
            from build_tf_gene_network import load_rna_counts, load_metadata, compute_rna_de
            rna_counts = load_rna_counts(rna_counts_path)
            metadata = load_metadata(metadata_path)
            contrast_cfg = cfg['contrasts'].get(contrast, cfg['contrasts']['NvT'])
            rna_de = compute_rna_de(rna_counts, metadata, contrast_cfg, cfg['scoring'])
        else:
            print("ERROR: Cannot find RNA data files")
            sys.exit(1)
    else:
        print(f"\nLoading RNA DE from: {rna_de_path}")
        rna_de = pd.read_csv(rna_de_path, sep='\t')

    print(f"  Loaded {len(rna_de)} genes with DE results")

    # Load peak2gene mapping
    peak2gene_path = base_dir / cfg['output']['base_dir'] / contrast / cfg['output']['tables_dir'] / f"peak2gene_{contrast}.tsv"

    if not peak2gene_path.exists():
        print(f"WARNING: Peak2gene mapping not found: {peak2gene_path}")
        print("  Attempting to create peak2gene mapping...")
        # Fall back to creating from raw data
        from build_tf_gene_network import load_peak_annotations, create_peak2gene_mapping
        peak_annot_path = resolve_path(cfg['data']['atac_peaks_annotated'], base_dir)
        peak_annot = load_peak_annotations(peak_annot_path)
        peak2gene = create_peak2gene_mapping(peak_annot, cfg)
    else:
        print(f"Loading peak2gene from: {peak2gene_path}")
        peak2gene = pd.read_csv(peak2gene_path, sep='\t')

    print(f"  Loaded {len(peak2gene)} peak-gene links")

    # Setup output directory
    out_dir = Path(args.output) if args.output else base_dir / cfg['output']['base_dir'] / "tf_focus"
    out_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = out_dir / cfg['output'].get('tables_dir', 'tables')
    tables_dir.mkdir(exist_ok=True)

    print(f"\nOutput directory: {out_dir}")

    # Process each TF
    for tf_info in tfs_to_process:
        tf_name = tf_info['name']
        motif_id = tf_info['motif_id']

        print(f"\n{'='*60}")
        print(f"Processing: {tf_name} ({motif_id})")
        print("="*60)

        # Build cascade network
        edges_df, nodes_df = build_cascade_network(
            primary_tf=tf_name,
            primary_motif=motif_id,
            tf_list=tf_list,
            peak2gene=peak2gene,
            rna_de=rna_de,
            bindetect_dir=bindetect_dir,
            config=tf_cfg
        )

        if len(edges_df) == 0:
            print(f"  WARNING: No edges generated for {tf_name}, skipping")
            continue

        # Save outputs
        contrast_label = tf_cfg.get('contrast', {}).get('name', 'contrast').replace(' ', '_')
        edges_file = tables_dir / f"tf_focus_edges_{tf_name}_{contrast_label}.tsv"
        nodes_file = tables_dir / f"tf_focus_nodes_{tf_name}_{contrast_label}.tsv"

        edges_df.to_csv(edges_file, sep='\t', index=False)
        nodes_df.to_csv(nodes_file, sep='\t', index=False)

        print(f"\n  Saved: {edges_file}")
        print(f"  Saved: {nodes_file}")

    print("\n" + "="*60)
    print("Done!")
    print("="*60)


if __name__ == '__main__':
    main()
