# Workflow Guide

This guide walks through the complete workflow for using RNA-ATAC-TF-Network.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Step 1: Prepare Input Data](#step-1-prepare-input-data)
3. [Step 2: Configure the Pipeline](#step-2-configure-the-pipeline)
4. [Step 3: Run the General Network Analysis](#step-3-run-the-general-network-analysis)
5. [Step 4: Run TF-Focused Analysis](#step-4-run-tf-focused-analysis)
6. [Step 5: Interpret Outputs](#step-5-interpret-outputs)
7. [Common Workflows](#common-workflows)
8. [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Required Software

- Python 3.8 or higher
- R 4.0 or higher
- TOBIAS (for generating input data)

### Installation

```bash
# Clone the repository
git clone https://github.com/your-username/RNA-ATAC-TF-Network.git
cd RNA-ATAC-TF-Network

# Run setup script (installs Python and R dependencies)
./install/setup.sh
```

### Required Input Data

You need outputs from:
1. **TOBIAS**: TF footprinting and BINDetect differential binding analysis
2. **ATAC-seq**: Peak calling and annotation (e.g., with ChIPseeker)
3. **RNA-seq**: Gene expression quantification (e.g., with featureCounts)

See [INPUT_DATA_FORMAT.md](INPUT_DATA_FORMAT.md) for detailed specifications.

---

## Step 1: Prepare Input Data

### From TOBIAS

Run TOBIAS on your ATAC-seq data:

```bash
# 1. Correct Tn5 bias
TOBIAS ATACorrect --bam sample.bam --genome genome.fa --peaks peaks.bed --outdir corrected/

# 2. Calculate footprint scores
TOBIAS ScoreBigwig --signal corrected/sample_corrected.bw --regions peaks.bed --output sample_footprints.bw

# 3. Run BINDetect for differential binding
TOBIAS BINDetect --motifs jaspar_motifs.meme --signals *.bw --genome genome.fa \
    --peaks peaks.bed --outdir BINDetect/ --cond_names Condition1 Condition2
```

This produces:
- `BINDetect/bindetect_results.txt` - Differential binding results
- `BINDetect/TF_NAME_MOTIF/beds/` - Per-TF binding site BED files

### From ATAC-seq

Annotate your peaks with ChIPseeker in R:

```r
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)  # or appropriate TxDb

peaks <- readPeakFile("peaks.bed")
peakAnno <- annotatePeak(peaks, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
write.table(as.data.frame(peakAnno), "peaks_annotated.txt", sep="\t", quote=FALSE)
```

### From RNA-seq

Generate count matrix with featureCounts:

```bash
featureCounts -a annotation.gtf -o counts.txt *.bam
```

Create metadata CSV manually:
```csv
Probe,Condition,CellType
Sample1,Normal,TypeA
Sample2,Normal,TypeA
Sample3,Transformed,TypeB
Sample4,Transformed,TypeB
```

---

## Step 2: Configure the Pipeline

### Create Configuration File

```bash
# Copy template
cp config/config_template.yaml config/config.yaml

# Edit with your settings
nano config/config.yaml
```

### Key Configuration Sections

#### Data Paths

```yaml
data:
  tf_binding_matrix: "/path/to/TF_MAT.csv"
  bindetect_results: "/path/to/bindetect_results.txt"
  atac_peaks_annotated: "/path/to/peaks_annotated.txt"
  atac_gene_counts: "/path/to/atac_gene_counts.txt"
  rna_counts: "/path/to/rna_counts.csv"
  metadata: "/path/to/metadata.csv"
```

#### Define Contrasts

```yaml
contrasts:
  TreatmentVsControl:
    name: "Treatment_vs_Control"
    description: "Compare treated to control samples"
    atac_groups:
      case: ["Treatment_rep1", "Treatment_rep2"]
      control: ["Control_rep1", "Control_rep2"]
    rna_filter:
      column: "Condition"
      case: "Treated"
      control: "Control"
```

#### Focus Pathways

```yaml
msigdb:
  species: "Mus musculus"  # or "Homo sapiens"
  focus_pathways:
    - "HALLMARK_INTERFERON_ALPHA_RESPONSE"
    - "HALLMARK_INFLAMMATORY_RESPONSE"
    # Add pathways relevant to your biology
```

#### TF-Focused Settings (Optional)

```yaml
tf_focus:
  enabled: true
  bindetect_bed_dir: "/path/to/BINDetect/"
  tfs:
    - name: "STAT1"
      motif_id: "MA0137.3"
    - name: "CEBPA"
      motif_id: "MA0102.5"
```

---

## Step 3: Run the General Network Analysis

### Basic Run

```bash
# Run for a single contrast
./run_pipeline.sh TreatmentVsControl
```

### What Happens

1. **Python script** (`build_tf_gene_network.py`):
   - Loads all input data
   - Creates peak-to-gene mapping
   - Computes TF activity scores
   - Computes accessibility changes
   - Computes RNA differential expression
   - Builds TF→gene edges with integrated scores
   - Outputs tables to `output/<contrast>/tables/`

2. **R script** (`plot_pathway_network.R`):
   - Loads MSigDB gene sets
   - Builds pathway-specific networks
   - Generates network visualizations
   - Creates enrichment maps
   - Exports Cytoscape files

### Run All Contrasts

```bash
./run_pipeline.sh all
```

### Check Progress

The pipeline prints progress messages:
```
============================================================
Building TF-Gene Network for: Treatment_vs_Control
============================================================

Loading BINDetect results from bindetect_results.txt...
    Loaded 600 TFs with 25 columns

Computing TF activity scores...
    Active TFs (|activity| >= 0.01): 245

Building TF→gene edges...
    Final edges: 15234 (245 TFs, 2341 genes)
```

---

## Step 4: Run TF-Focused Analysis

### Enable in Config

```yaml
tf_focus:
  enabled: true
  tfs:
    - name: "CEBPA"
      motif_id: "MA0102.5"
```

### Run TF-Focused Analysis

```bash
# Run TF-focus only
./run_pipeline.sh tf_focus

# Run with main pipeline
./run_pipeline.sh TreatmentVsControl --tf_focus
```

### What Happens

1. **Python script** (`build_tf_focus_network.py`):
   - Loads TF-specific binding sites from BED files
   - Maps binding sites to target genes
   - Identifies secondary TFs among targets
   - Builds cascade network (Layer 0 → 1 → 2)
   - Outputs tables to `output/tf_focus/tables/`

2. **R script** (`plot_tf_focus_network.R`):
   - Creates hierarchical cascade visualizations
   - Generates circular layout versions
   - Creates compact versions for figure panels

---

## Step 5: Interpret Outputs

### Output Directory Structure

```
output/
├── TreatmentVsControl/
│   ├── tables/
│   │   ├── tf_gene_edges_TreatmentVsControl.tsv
│   │   ├── tf_gene_nodes_TreatmentVsControl.tsv
│   │   ├── tf_activity_TreatmentVsControl.tsv
│   │   ├── rna_de_TreatmentVsControl.tsv
│   │   └── peak2gene_TreatmentVsControl.tsv
│   ├── plots/
│   │   ├── network_TreatmentVsControl_HALLMARK_*.png
│   │   ├── network_TreatmentVsControl_HALLMARK_*.pdf
│   │   └── emap_TreatmentVsControl.png
│   └── cytoscape/
│       ├── cytoscape_edges_TreatmentVsControl.tsv
│       └── cytoscape_nodes_TreatmentVsControl.tsv
└── tf_focus/
    ├── tables/
    │   ├── tf_focus_edges_CEBPA_*.tsv
    │   └── tf_focus_nodes_CEBPA_*.tsv
    └── plots/
        ├── tf_focus_CEBPA.png
        ├── tf_focus_CEBPA_circular.png
        └── tf_focus_CEBPA_compact.png
```

### Key Output Files

#### Edge Table (`tf_gene_edges_*.tsv`)

| Column | Description |
|--------|-------------|
| `tf` | Transcription factor name |
| `gene` | Target gene name |
| `edge_score` | Integrated regulatory score |
| `tf_activity` | TF differential binding score |
| `atac_log2fc` | Accessibility change |
| `rna_log2fc` | Expression change |
| `rna_pvalue` | Expression significance |

#### Node Table (`tf_gene_nodes_*.tsv`)

| Column | Description |
|--------|-------------|
| `node` | Node name (TF or gene) |
| `node_type` | "TF" or "gene" |
| `tf_activity` | TF activity (for TFs) |
| `expression` | log2FC (for genes) |
| `degree` | Number of connections |

### Visual Interpretation

**Network Plots**:
- Diamond shapes = TFs
- Circles = genes
- Node color = expression/activity (red=up, blue=down)
- Edge width = regulatory strength

**Cascade Plots**:
- Center = primary TF
- Inner ring = direct targets (Layer 1)
- Outer ring = secondary targets (Layer 2)

---

## Common Workflows

### Workflow 1: Identify Key Regulatory TFs

1. Run general network analysis
2. Check `tf_activity_*.tsv` for highly active TFs
3. Examine pathway-specific plots
4. Run TF-focused analysis on top candidates

### Workflow 2: Compare Two Conditions

1. Define contrasts in config
2. Run both contrasts:
   ```bash
   ./run_pipeline.sh Condition1
   ./run_pipeline.sh Condition2
   ```
3. Compare edge tables and TF activities

### Workflow 3: Explore in Cytoscape

1. Run pipeline to generate Cytoscape files
2. Import `cytoscape_edges_*.tsv` as network
3. Import `cytoscape_nodes_*.tsv` as node attributes
4. Apply visual styles based on node_type, expression

---

## Troubleshooting

### "No edges generated"

**Cause**: Filters too stringent or data issues.

**Solutions**:
- Lower `rna_log2fc_min` threshold in config
- Lower `tf_change_min` threshold
- Check that gene names match between files
- Verify coordinate systems (0-based vs 1-based)

### "No TFs found for pathway"

**Cause**: No TFs regulate genes in the pathway.

**Solutions**:
- Try broader pathways (Hallmark instead of KEGG)
- Lower `min_tf_activity` threshold
- Check that TF names in bindetect match expected format

### "BINDetect directory not found"

**Cause**: Path to TOBIAS BINDetect output is incorrect.

**Solutions**:
- Verify `bindetect_bed_dir` path in config
- Check that directory contains TF_MOTIF subdirectories
- Ensure BED files exist in `*/beds/` subdirectories

### "Memory error"

**Cause**: Large datasets exceeding available RAM.

**Solutions**:
- Increase `rna_log2fc_min` to analyze fewer genes
- Reduce `max_genes_per_pathway`
- Process contrasts one at a time
- Use a machine with more RAM

### Plots look wrong

**Cause**: Too many or too few nodes.

**Solutions**:
- Adjust `max_genes_per_pathway` and `max_tfs_per_pathway`
- Adjust `edge_score_min` to filter weak edges
- Try different layout algorithms
