# Input Data Format Specifications

This document describes the required format for each input file used by RNA-ATAC-TF-Network.

## Table of Contents

1. [TOBIAS Outputs](#tobias-outputs)
   - [TF Binding Matrix](#tf-binding-matrix-tf_matcsv)
   - [BINDetect Results](#bindetect-results-bindetect_resultstxt)
   - [BINDetect BED Files](#bindetect-bed-files)
2. [ATAC-seq Data](#atac-seq-data)
   - [Peak Annotations](#peak-annotations-peaks_annotatedtxt)
   - [Gene-level ATAC Counts](#gene-level-atac-counts)
   - [Peak-level ATAC Counts](#peak-level-atac-counts)
3. [RNA-seq Data](#rna-seq-data)
   - [Count Matrix](#rna-count-matrix-rna_countscsv)
   - [Metadata](#sample-metadata-metadatacsv)

---

## TOBIAS Outputs

### TF Binding Matrix (`TF_MAT.csv`)

A matrix of TF binding counts/scores across samples, output from TOBIAS.

**Format**: CSV with TF names as row names and sample names as columns.

**Example**:
```csv
,Sample1_ctrl,Sample1_treat,Sample2_ctrl,Sample2_treat
STAT1,1523,2847,1612,2953
STAT3,982,1834,1045,1892
IRF1,756,1423,812,1567
CEBPA,2341,1876,2456,1923
```

**Requirements**:
- First column (unnamed or "TF") contains TF names
- Column names match your sample identifiers
- Values are numeric (counts or scores)

---

### BINDetect Results (`bindetect_results.txt`)

Differential TF binding results from TOBIAS BINDetect.

**Format**: Tab-separated file with one row per TF.

**Required columns**:
- `name`: TF name
- `*_mean_score` columns: Mean footprint scores per condition
- `*_change` columns: Pairwise differential binding scores
- `*_pvalue` columns: P-values for differential binding

**Example**:
```
name	motif_id	"Condition1_ctrl"_mean_score	"Condition1_treat"_mean_score	"Condition1_ctrl"_"Condition1_treat"_change	"Condition1_ctrl"_"Condition1_treat"_pvalue
STAT1	MA0137.3	0.234	0.456	0.222	0.001
STAT3	MA0144.2	0.187	0.298	0.111	0.023
```

**Notes**:
- Column names may contain Unicode smart quotes (`"..."`) or regular quotes
- The pipeline handles both formats automatically

---

### BINDetect BED Files

Directory structure containing per-TF binding site BED files.

**Structure**:
```
BINDetect/
├── STAT1_MA0137.3/
│   └── beds/
│       ├── STAT1_MA0137.3_all.bed
│       ├── STAT1_MA0137.3_"Condition1_ctrl"_bound.bed
│       ├── STAT1_MA0137.3_"Condition1_treat"_bound.bed
│       └── ...
├── STAT3_MA0144.2/
│   └── beds/
│       └── ...
```

**BED file format** (17 columns):
| Column | Description |
|--------|-------------|
| 1-3 | Motif site coordinates (chr, start, end) |
| 4 | TF name |
| 5 | Binding score (motif affinity) |
| 6 | Strand |
| 7-9 | Peak coordinates (chr, start, end) |
| 10-17 | Accessibility values for each condition |

**Example**:
```
chr1	1000	1010	STAT1	8.5	+	chr1	900	1100	0.23	0.45	0.21	0.48	0.19	0.42	0.25	0.51
chr2	5000	5010	STAT1	7.2	-	chr2	4800	5200	0.18	0.39	0.16	0.41	0.15	0.38	0.19	0.44
```

---

## ATAC-seq Data

### Peak Annotations (`peaks_annotated.txt`)

Peak annotations with gene assignments, typically from ChIPseeker or HOMER.

**Format**: Tab-separated file.

**Required columns**:
| Column | Description |
|--------|-------------|
| `seqnames` | Chromosome |
| `start` | Peak start position |
| `end` | Peak end position |
| `geneId` | Assigned gene symbol |
| `distanceToTSS` | Distance to transcription start site |
| `annotation` | Annotation type (Promoter, Intron, Distal, etc.) |

**Example**:
```
seqnames	start	end	width	strand	geneId	distanceToTSS	annotation
chr1	1000	1500	501	*	Gene1	-250	Promoter (<=1kb)
chr1	5000	5800	801	*	Gene2	25000	Distal Intergenic
chr2	3000	3600	601	*	Gene3	500	Promoter (<=1kb)
chr2	8000	8400	401	*	Gene4	-5000	Intron (ENSG00000...)
```

**Notes**:
- The `annotation` column is used to classify peaks as promoter, intronic, or distal
- Promoter peaks get higher weights in the analysis

---

### Gene-level ATAC Counts

Gene-level accessibility counts, summarized from peak counts.

**Format**: Tab-separated file with gene names as row names.

**Example**:
```
Gene	Sample1_ctrl	Sample1_treat	Sample2_ctrl	Sample2_treat
Gene1	245	412	267	398
Gene2	89	156	102	148
Gene3	523	489	545	472
```

---

### Peak-level ATAC Counts

Raw read counts per peak, used for more detailed analysis.

**Format**: Tab-separated file with peak IDs as row names.

**Example**:
```
peak_id	Sample1_ctrl	Sample1_treat	Sample2_ctrl	Sample2_treat
chr1:1000-1500	123	256	134	248
chr1:5000-5800	67	89	72	94
chr2:3000-3600	234	198	245	187
```

---

## RNA-seq Data

### RNA Count Matrix (`rna_counts.csv`)

Gene expression count matrix from featureCounts, Salmon, or similar.

**Format**: CSV with gene symbols as row names.

**Example**:
```csv
,Sample1_ctrl,Sample1_treat,Sample2_ctrl,Sample2_treat
Gene1,1523,2847,1612,2953
Gene2,982,1834,1045,1892
Gene3,756,1423,812,1567
```

**Requirements**:
- First column contains gene symbols
- Column names match your sample identifiers
- Values are raw counts (not normalized)

**Optional columns**:
- `Description`: Gene description (will be ignored)

---

### Sample Metadata (`metadata.csv`)

Sample information including conditions and groups.

**Format**: CSV file.

**Required columns**:
| Column | Description |
|--------|-------------|
| `Probe` or `Sample` | Sample identifier (matches RNA column names) |
| `Condition` | Primary condition grouping |

**Optional columns**:
| Column | Description |
|--------|-------------|
| `CellType` | Cell type or line |
| `Stimulus` | Treatment/stimulus |
| `Batch` | Batch information |

**Example**:
```csv
Probe,Condition,CellType,Stimulus
Sample1_ctrl,Normal,TypeA,Control
Sample1_treat,Normal,TypeA,Treated
Sample2_ctrl,Transformed,TypeB,Control
Sample2_treat,Transformed,TypeB,Treated
```

**Notes**:
- The `Condition` column is used for defining contrasts
- Sample names should match (with optional prefix differences) between RNA counts and metadata
- The pipeline handles common prefixes like 'r' (e.g., `rSample1` in counts maps to `Sample1` in metadata)

---

## Tips for Data Preparation

### Coordinate Systems

Be aware of coordinate systems:
- **BED format**: 0-based, half-open intervals `[start, end)`
- **GTF/GFF format**: 1-based, closed intervals `[start, end]`

The pipeline attempts to match coordinates automatically, but inconsistencies can cause mapping failures.

### Gene Naming

Use consistent gene naming:
- Prefer official gene symbols (e.g., `STAT1`, not `ENSG00000115415`)
- Ensure gene names match between RNA counts and peak annotations
- Case may matter - standardize to UPPERCASE or mixed case consistently

### Sample Naming

Use clear, consistent sample names:
- Include condition and replicate information
- Avoid special characters (except underscores)
- Good: `Condition1_ctrl_rep1`, `Treatment_IL10_3`
- Bad: `Sample #1`, `ctrl(new)`

### Quality Control

Before running the pipeline:
1. Verify file formats match specifications
2. Check for missing values or empty cells
3. Confirm sample names align across all files
4. Validate that gene symbols are consistent
