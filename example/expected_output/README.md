# Expected Output Structure

After running the pipeline, you should see the following output structure:

```
output/
├── <contrast_name>/
│   ├── tables/
│   │   ├── tf_gene_edges_<contrast>.tsv      # TF→gene edges with scores
│   │   ├── tf_gene_nodes_<contrast>.tsv      # Node attributes
│   │   ├── tf_activity_<contrast>.tsv        # TF activity scores
│   │   ├── rna_de_<contrast>.tsv             # RNA differential expression
│   │   └── peak2gene_<contrast>.tsv          # Peak-gene mapping
│   │
│   ├── plots/
│   │   ├── network_<contrast>_HALLMARK_*.png # Pathway network plots
│   │   ├── network_<contrast>_HALLMARK_*.pdf
│   │   ├── emap_<contrast>.png               # Enrichment map
│   │   └── emap_<contrast>.pdf
│   │
│   └── cytoscape/
│       ├── cytoscape_edges_<contrast>.tsv    # For Cytoscape import
│       └── cytoscape_nodes_<contrast>.tsv
│
└── tf_focus/
    ├── tables/
    │   ├── tf_focus_edges_<TF>_<contrast>.tsv
    │   └── tf_focus_nodes_<TF>_<contrast>.tsv
    │
    └── plots_<TF>_<contrast>/
        ├── tf_focus_<TF>.png                 # Hierarchical layout
        ├── tf_focus_<TF>.pdf
        ├── tf_focus_<TF>_circular.png        # Radial layout
        ├── tf_focus_<TF>_circular.pdf
        ├── tf_focus_<TF>_compact.png         # For figure panels
        └── tf_focus_<TF>_compact.pdf
```

## Example Files

### Edge Table Example (`tf_gene_edges_*.tsv`)

```
tf	gene	edge_score	edge_score_norm	n_peaks	tf_activity	atac_log2fc	rna_log2fc	rna_pvalue
STAT1	Irf1	0.542	0.876	3	0.234	0.89	2.34	0.001
STAT1	Gbp4	0.489	0.791	2	0.234	0.67	1.98	0.003
IRF1	Cxcl10	0.412	0.666	1	0.187	0.45	3.21	0.0001
```

### Node Table Example (`tf_gene_nodes_*.tsv`)

```
node	node_type	tf_activity	expression	direction	degree
STAT1	TF	0.234	NA	1	45
IRF1	TF	0.187	NA	1	32
Irf1	gene	NA	2.34	1	3
Gbp4	gene	NA	1.98	1	2
```

### TF Activity Example (`tf_activity_*.tsv`)

```
tf	tf_activity	tf_pvalue	direction
STAT1	0.234	0.001	1
STAT3	0.156	0.012	1
IRF1	0.187	0.003	1
NFKB1	-0.089	0.045	-1
```

## Verification Checklist

After running the pipeline, verify:

- [ ] Edge file contains expected TFs and genes
- [ ] Node file has both TF and gene entries
- [ ] TF activity values are reasonable (typically -1 to 1)
- [ ] RNA log2FC values match expected direction of change
- [ ] Network plots are generated for focus pathways
- [ ] Cytoscape files can be imported successfully
- [ ] TF cascade plots show correct hierarchy (if tf_focus enabled)

## Troubleshooting Empty Outputs

If output files are empty or missing:

1. Check log messages for errors
2. Verify input file paths in config.yaml
3. Lower scoring thresholds (rna_log2fc_min, tf_change_min)
4. Verify gene names match between RNA and ATAC data
5. Check that contrast group names match between config and data
