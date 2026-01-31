#!/usr/bin/env Rscript
################################################################################
#  cnet-style TF-Gene Regulatory Network Plots
#  ============================================
#  Creates publication-quality cnet-style plots similar to enrichplot::cnetplot
#  but designed for TF→gene regulatory networks.
#
#  This script focuses on:
#  - TF hubs with gene targets radiating outward
#  - Color by expression/activity
#  - Edge weights by regulatory score
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(RColorBrewer)
  library(scales)
  library(stringr)
  library(tidygraph)
})

# ==============================================================================
# CNET-STYLE PLOTTING FUNCTION
# ==============================================================================

plot_cnet_tf_gene <- function(edges,
                               nodes,
                               pathway_genes = NULL,
                               pathway_name = NULL,
                               top_n_tfs = 8,
                               top_n_genes_per_tf = 10,
                               min_edge_score = 0.1,
                               show_all_genes = FALSE,
                               tf_label_size = 4,
                               gene_label_size = 2.5,
                               layout = "kk",  # "kk", "fr", "star", "circle"
                               color_by = "expression",  # "expression" or "edge_score"
                               title = NULL) {
  #' Create cnet-style TF-gene regulatory network plot
  #'
  #' @param edges Edge table with columns: tf, gene, edge_score, rna_log2fc
  #' @param nodes Node table with columns: node, node_type, tf_activity, expression
  #' @param pathway_genes Optional: subset to genes in this pathway
  #' @param pathway_name Optional: name for title
  #' @param top_n_tfs Number of TFs to show (by number of targets)
  #' @param top_n_genes_per_tf Max genes to show per TF
  #' @param min_edge_score Minimum edge score threshold
  #' @param show_all_genes If TRUE, show shared genes between TFs
  #' @param layout ggraph layout algorithm
  #' @param color_by Color nodes by "expression" or "edge_score"

  # Filter to pathway if specified
  if (!is.null(pathway_genes)) {
    edges <- edges %>% filter(gene %in% pathway_genes)
  }

  if (nrow(edges) == 0) {
    message("No edges to plot")
    return(NULL)
  }

  # Filter by edge score
  edges <- edges %>% filter(abs(edge_score) >= min_edge_score)

  # Get top TFs by number of targets
  tf_stats <- edges %>%
    group_by(tf) %>%
    summarise(
      n_targets = n(),
      mean_score = mean(abs(edge_score)),
      .groups = "drop"
    ) %>%
    arrange(desc(n_targets))

  top_tfs <- head(tf_stats$tf, top_n_tfs)

  # Filter edges to top TFs
  edges_plot <- edges %>%
    filter(tf %in% top_tfs)

  # For each TF, keep top genes by edge score
  if (!show_all_genes) {
    edges_plot <- edges_plot %>%
      group_by(tf) %>%
      slice_max(order_by = abs(edge_score), n = top_n_genes_per_tf, with_ties = FALSE) %>%
      ungroup()
  }

  # Get unique nodes
  all_tfs <- unique(edges_plot$tf)
  all_genes <- unique(edges_plot$gene)

  # Build node data frame
  node_df <- bind_rows(
    # TF nodes
    data.frame(
      node = all_tfs,
      type = "TF",
      stringsAsFactors = FALSE
    ) %>%
      left_join(
        nodes %>% filter(node_type == "TF") %>% select(node, tf_activity),
        by = "node"
      ) %>%
      mutate(value = tf_activity),

    # Gene nodes
    data.frame(
      node = all_genes,
      type = "gene",
      stringsAsFactors = FALSE
    ) %>%
      left_join(
        edges_plot %>%
          group_by(gene) %>%
          summarise(rna_log2fc = first(rna_log2fc), .groups = "drop") %>%
          rename(node = gene),
        by = "node"
      ) %>%
      mutate(value = rna_log2fc)
  )

  # Create igraph
  g <- graph_from_data_frame(
    d = edges_plot %>% select(tf, gene, edge_score, rna_log2fc),
    directed = TRUE,
    vertices = node_df
  )

  # Convert to tidygraph for easier manipulation
  tg <- as_tbl_graph(g)

  # Calculate node degrees
  tg <- tg %>%
    activate(nodes) %>%
    mutate(
      degree = centrality_degree(mode = "all"),
      is_tf = type == "TF"
    )

  # Determine layout
  if (layout == "star") {
    # Star layout with TFs in center
    layout_fn <- create_layout(tg, layout = "star", center = which(V(g)$type == "TF")[1])
  } else {
    layout_fn <- create_layout(tg, layout = layout)
  }

  # Color scale limits
  val_range <- range(node_df$value, na.rm = TRUE)
  val_max <- max(abs(val_range))

  # Create plot
  p <- ggraph(layout_fn) +
    # Edges
    geom_edge_link(
      aes(
        edge_alpha = abs(edge_score),
        edge_width = abs(edge_score)
      ),
      color = "grey50",
      arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
      end_cap = circle(2, "mm")
    ) +
    scale_edge_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    scale_edge_width_continuous(range = c(0.3, 1.5), guide = "none") +

    # Gene nodes (circles)
    geom_node_point(
      data = . %>% filter(!is_tf),
      aes(size = degree, fill = value),
      shape = 21, color = "grey30", stroke = 0.3
    ) +

    # TF nodes (larger, diamonds)
    geom_node_point(
      data = . %>% filter(is_tf),
      aes(size = degree * 2, fill = value),
      shape = 23, color = "black", stroke = 1
    ) +

    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-val_max, val_max),
      name = ifelse(color_by == "expression", "log2FC", "Score"),
      na.value = "grey70"
    ) +
    scale_size_continuous(range = c(2, 12), guide = "none") +

    # TF labels (bold, larger)
    geom_node_text(
      data = . %>% filter(is_tf),
      aes(label = name),
      size = tf_label_size, fontface = "bold",
      repel = TRUE,
      point.padding = unit(0.5, "lines"),
      box.padding = unit(0.3, "lines"),
      max.overlaps = 20
    ) +

    # Gene labels (smaller)
    geom_node_text(
      data = . %>% filter(!is_tf),
      aes(label = name),
      size = gene_label_size,
      repel = TRUE,
      point.padding = unit(0.2, "lines"),
      box.padding = unit(0.2, "lines"),
      max.overlaps = 30,
      segment.color = "grey70",
      segment.size = 0.2
    ) +

    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    )

  # Add title
  if (!is.null(title)) {
    p <- p + labs(title = str_wrap(title, 60))
  } else if (!is.null(pathway_name)) {
    clean_name <- gsub("_", " ", gsub("^HALLMARK_|^KEGG_|^REACTOME_", "", pathway_name))
    p <- p + labs(title = str_wrap(clean_name, 60))
  }

  # Add subtitle with stats
  p <- p + labs(
    subtitle = sprintf("%d TFs regulating %d genes",
                       length(all_tfs), length(all_genes))
  )

  return(p)
}

# ==============================================================================
# MULTI-PATHWAY COMPARISON PLOT
# ==============================================================================

plot_tf_pathway_heatmap <- function(edges, nodes, gene_sets,
                                     pathways = NULL,
                                     top_n_pathways = 15,
                                     top_n_tfs = 20,
                                     min_genes = 3) {
  #' Create heatmap of TF activity across pathways
  #'
  #' Shows which TFs regulate genes in which pathways

  # Calculate TF activity per pathway
  regulated_genes <- unique(edges$gene)

  if (is.null(pathways)) {
    # Find pathways with most regulated genes
    pw_stats <- sapply(names(gene_sets), function(pw) {
      length(intersect(gene_sets[[pw]], regulated_genes))
    })
    pathways <- names(sort(pw_stats, decreasing = TRUE))[1:top_n_pathways]
    pathways <- pathways[pw_stats[pathways] >= min_genes]
  }

  # For each TF, count regulated genes per pathway
  tfs <- unique(edges$tf)

  tf_pathway_mat <- sapply(pathways, function(pw) {
    pw_genes <- gene_sets[[pw]]
    sapply(tfs, function(tf) {
      tf_genes <- edges$gene[edges$tf == tf]
      length(intersect(tf_genes, pw_genes))
    })
  })

  # Filter to TFs with at least some pathway hits
  tf_sums <- rowSums(tf_pathway_mat)
  keep_tfs <- names(sort(tf_sums, decreasing = TRUE))[1:min(top_n_tfs, sum(tf_sums > 0))]
  tf_pathway_mat <- tf_pathway_mat[keep_tfs, , drop = FALSE]

  # Convert to long format for ggplot
  mat_df <- as.data.frame(tf_pathway_mat) %>%
    rownames_to_column("TF") %>%
    pivot_longer(-TF, names_to = "Pathway", values_to = "n_genes")

  # Clean pathway names
  mat_df$Pathway_clean <- sapply(mat_df$Pathway, function(x) {
    x <- gsub("^HALLMARK_|^KEGG_|^REACTOME_|^GO_", "", x)
    x <- gsub("_", " ", x)
    str_wrap(x, 30)
  })

  # Get TF activity for color
  tf_act <- nodes %>%
    filter(node_type == "TF") %>%
    select(node, tf_activity) %>%
    rename(TF = node)

  mat_df <- mat_df %>%
    left_join(tf_act, by = "TF")

  # Order TFs by activity
  tf_order <- tf_act %>%
    filter(TF %in% keep_tfs) %>%
    arrange(desc(tf_activity)) %>%
    pull(TF)

  mat_df$TF <- factor(mat_df$TF, levels = tf_order)

  # Order pathways by total regulation
  pw_order <- mat_df %>%
    group_by(Pathway_clean) %>%
    summarise(total = sum(n_genes), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(Pathway_clean)

  mat_df$Pathway_clean <- factor(mat_df$Pathway_clean, levels = rev(pw_order))

  # Create heatmap
  p <- ggplot(mat_df, aes(x = TF, y = Pathway_clean, fill = n_genes)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(
      low = "white", high = "#B2182B",
      name = "Regulated\nGenes"
    ) +
    geom_text(
      aes(label = ifelse(n_genes > 0, n_genes, "")),
      size = 2.5, color = "grey20"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = "TF Regulation Across Pathways")

  return(p)
}

# ==============================================================================
# EXAMPLE USAGE / DEMO
# ==============================================================================

demo_plot <- function() {
  # This function demonstrates the plotting with example data
  # Useful for testing without running the full pipeline

  message("Creating demo plot with synthetic data...")

  # Create synthetic edges
  tfs <- c("STAT1", "STAT3", "IRF1", "IRF7", "NFKB1", "JUN", "FOS", "MYC")
  genes_per_tf <- list(
    STAT1 = c("OAS1", "OAS2", "MX1", "MX2", "IFIT1", "IFIT2", "ISG15", "IRF9"),
    STAT3 = c("BCL2", "MCL1", "SOCS3", "CCND1", "MYC", "VEGFA", "HIF1A"),
    IRF1 = c("STAT1", "IRF7", "CXCL10", "CXCL9", "GBP1", "GBP2"),
    IRF7 = c("IFNA1", "IFNB1", "IRF3", "OAS1", "ISG15"),
    NFKB1 = c("IL6", "TNF", "IL1B", "CXCL8", "CCL2", "ICAM1"),
    JUN = c("FOS", "JUNB", "ATF3", "CDKN1A", "GADD45A"),
    FOS = c("JUN", "FOSB", "EGR1", "MMP1"),
    MYC = c("CDK4", "CCND1", "LDHA", "PKM", "ODC1")
  )

  edges <- do.call(rbind, lapply(names(genes_per_tf), function(tf) {
    genes <- genes_per_tf[[tf]]
    data.frame(
      tf = tf,
      gene = genes,
      edge_score = runif(length(genes), 0.2, 1) * sample(c(-1, 1), length(genes), replace = TRUE),
      rna_log2fc = rnorm(length(genes), mean = 0.5, sd = 1),
      stringsAsFactors = FALSE
    )
  }))

  # Create nodes
  all_tfs <- unique(edges$tf)
  all_genes <- unique(edges$gene)

  nodes <- rbind(
    data.frame(
      node = all_tfs,
      node_type = "TF",
      tf_activity = rnorm(length(all_tfs), 0.5, 0.3),
      expression = NA,
      stringsAsFactors = FALSE
    ),
    data.frame(
      node = all_genes,
      node_type = "gene",
      tf_activity = NA,
      expression = edges$rna_log2fc[match(all_genes, edges$gene)],
      stringsAsFactors = FALSE
    )
  )

  # Create plot
  p <- plot_cnet_tf_gene(
    edges = edges,
    nodes = nodes,
    top_n_tfs = 6,
    top_n_genes_per_tf = 8,
    title = "Demo: TF-Gene Regulatory Network"
  )

  print(p)

  return(p)
}

# Run demo if called directly
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) > 0 && args[1] == "--demo") {
    p <- demo_plot()
    ggsave("demo_cnet_plot.pdf", p, width = 12, height = 10)
    message("Saved: demo_cnet_plot.pdf")
  } else {
    message("Usage: Rscript plot_cnet_style.R --demo")
    message("Or source this file in R and call the functions directly.")
  }
}
