#!/usr/bin/env Rscript
################################################################################
#  TF-Gene Regulatory Network Visualization
#  =========================================
#  Creates cnet/emap-style network plots for TF→gene regulatory relationships
#  integrated with MSigDB pathways.
#
#  Usage:
#    Rscript plot_pathway_network.R --edges edges.tsv --nodes nodes.tsv \
#      --contrast NvT --pathway HALLMARK_INTERFERON_ALPHA_RESPONSE
#
#    Rscript plot_pathway_network.R --edges edges.tsv --nodes nodes.tsv \
#      --contrast NvT --all_pathways
################################################################################

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(RColorBrewer)
  library(msigdbr)
  library(scales)
  library(stringr)
  library(viridis)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Default settings (can be overridden by config.yaml)
DEFAULT_CONFIG <- list(
  msigdb = list(
    species = "Mus musculus",
    collections = list(
      list(category = "H", subcategory = NULL),
      list(category = "C2", subcategory = "CP:KEGG"),
      list(category = "C2", subcategory = "CP:REACTOME")
    )
  ),
  network = list(
    max_genes_per_pathway = 50,
    max_tfs_per_pathway = 15,
    top_edges_per_tf = 20,
    min_tf_activity = 0.05
  ),
  output = list(
    plot_width = 12,
    plot_height = 10,
    plot_dpi = 300,
    save_pdf = TRUE,
    save_png = TRUE
  )
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

load_config <- function(config_path) {
  if (file.exists(config_path)) {
    cfg <- yaml::read_yaml(config_path)
    # Merge with defaults
    for (section in names(DEFAULT_CONFIG)) {
      if (!(section %in% names(cfg))) {
        cfg[[section]] <- DEFAULT_CONFIG[[section]]
      }
    }
    return(cfg)
  }
  return(DEFAULT_CONFIG)
}

wrap_label <- function(x, width = 40) {
  str_wrap(x, width = width)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_-]", "_", x)
}

# Load MSigDB gene sets
load_msigdb <- function(species = "Mus musculus", collections = NULL) {
  message("Loading MSigDB gene sets...")

  all_sets <- list()

  # Default collections if not specified
  if (is.null(collections)) {
    collections <- list(
      list(category = "H", subcategory = NULL),
      list(category = "C2", subcategory = "CP:KEGG")
    )
  }

  # Check msigdbr version and use appropriate API
  msigdbr_version <- packageVersion("msigdbr")
  use_new_api <- msigdbr_version >= "10.0.0"

  for (col in collections) {
    cat <- col$category
    subcat <- col$subcategory

    tryCatch({
      if (use_new_api) {
        # New API (msigdbr >= 10.0.0) uses collection/subcollection
        if (is.null(subcat) || subcat == "null") {
          sets <- msigdbr(species = species, collection = cat)
        } else {
          sets <- msigdbr(species = species, collection = cat, subcollection = subcat)
        }
      } else {
        # Old API uses category/subcategory
        if (is.null(subcat) || subcat == "null") {
          sets <- msigdbr(species = species, category = cat)
        } else {
          sets <- msigdbr(species = species, category = cat, subcategory = subcat)
        }
      }

      if (nrow(sets) > 0) {
        all_sets[[length(all_sets) + 1]] <- sets
        message(sprintf("  Loaded %s %s: %d gene sets",
                        cat, ifelse(is.null(subcat) || subcat == "null", "", subcat),
                        n_distinct(sets$gs_name)))
      }
    }, error = function(e) {
      message(sprintf("  Warning: Could not load %s %s: %s", cat,
                      ifelse(is.null(subcat), "", subcat), e$message))
    })
  }

  if (length(all_sets) == 0) {
    # Fall back to just Hallmark
    message("  Falling back to Hallmark only...")
    if (use_new_api) {
      sets <- msigdbr(species = species, collection = "H")
    } else {
      sets <- msigdbr(species = species, category = "H")
    }
    all_sets[[1]] <- sets
  }

  msig <- bind_rows(all_sets) %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  # Convert to list format
  gene_sets <- split(msig$gene_symbol, msig$gs_name)

  message(sprintf("  Total: %d gene sets", length(gene_sets)))
  return(gene_sets)
}

# ==============================================================================
# NETWORK BUILDING FUNCTIONS
# ==============================================================================

build_pathway_network <- function(edges, nodes, pathway_genes, pathway_name,
                                  max_genes = 50, max_tfs = 15,
                                  min_tf_activity = 0.05) {
  #' Build a network for a specific pathway
  #'
  #' @param edges TF-gene edge table

#' @param nodes Node attributes table
  #' @param pathway_genes Character vector of genes in the pathway
  #' @param pathway_name Name of the pathway
  #' @param max_genes Maximum number of genes to include
  #' @param max_tfs Maximum number of TFs to include

  # Filter edges to pathway genes
  pathway_edges <- edges %>%
    filter(gene %in% pathway_genes)

  if (nrow(pathway_edges) == 0) {
    message(sprintf("  No edges found for pathway: %s", pathway_name))
    return(NULL)
  }

  # Get TFs that regulate pathway genes
  pathway_tfs <- pathway_edges %>%
    group_by(tf) %>%
    summarise(
      n_targets = n(),
      mean_edge_score = mean(edge_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_targets))

  # Filter TFs by activity threshold
  tf_nodes <- nodes %>% filter(node_type == "TF")
  active_tfs <- tf_nodes %>%
    filter(abs(tf_activity) >= min_tf_activity) %>%
    pull(node)

  pathway_tfs <- pathway_tfs %>%
    filter(tf %in% active_tfs)

  if (nrow(pathway_tfs) == 0) {
    message(sprintf("  No active TFs found for pathway: %s", pathway_name))
    return(NULL)
  }

  # Limit TFs
  top_tfs <- head(pathway_tfs$tf, max_tfs)

  # Filter edges to top TFs
  pathway_edges <- pathway_edges %>%
    filter(tf %in% top_tfs)

  # Limit genes (by edge score)
  top_genes <- pathway_edges %>%
    group_by(gene) %>%
    summarise(max_score = max(abs(edge_score), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_score)) %>%
    head(max_genes) %>%
    pull(gene)

  pathway_edges <- pathway_edges %>%
    filter(gene %in% top_genes)

  # Build igraph - ensure numeric columns are proper vectors
  edge_df <- pathway_edges %>%
    select(tf, gene, edge_score, rna_log2fc) %>%
    mutate(
      edge_score = as.numeric(edge_score),
      rna_log2fc = as.numeric(rna_log2fc)
    )

  g <- graph_from_data_frame(
    d = edge_df,
    directed = TRUE,
    vertices = NULL
  )

  # Add node attributes
  V(g)$name <- V(g)$name
  V(g)$type <- ifelse(V(g)$name %in% top_tfs, "TF", "gene")

  # Add TF activity using a named vector lookup (more reliable)
  tf_act <- nodes %>%
    filter(node_type == "TF") %>%
    select(node, tf_activity) %>%
    distinct(node, .keep_all = TRUE)  # ensure unique nodes
  tf_act_vec <- setNames(tf_act$tf_activity, tf_act$node)
  V(g)$tf_activity <- ifelse(V(g)$name %in% names(tf_act_vec),
                              tf_act_vec[V(g)$name],
                              NA_real_)

  # Add gene expression using a named vector lookup
  gene_expr <- nodes %>%
    filter(node_type == "gene") %>%
    select(node, expression) %>%
    distinct(node, .keep_all = TRUE)
  gene_expr_vec <- setNames(gene_expr$expression, gene_expr$node)
  V(g)$expression <- ifelse(V(g)$name %in% names(gene_expr_vec),
                             gene_expr_vec[V(g)$name],
                             NA_real_)

  # Add degree
  V(g)$degree <- degree(g, mode = "all")

  list(
    graph = g,
    pathway_name = pathway_name,
    n_tfs = sum(V(g)$type == "TF"),
    n_genes = sum(V(g)$type == "gene"),
    n_edges = ecount(g)
  )
}

# ==============================================================================
# PLOTTING FUNCTIONS
# ==============================================================================

plot_tf_gene_network <- function(net, title = NULL,
                                  layout = "fr", # "fr", "kk", "stress", "circle"
                                  tf_size_range = c(6, 20),
                                  gene_size_range = c(3, 10),
                                  edge_alpha_range = c(0.3, 0.9)) {
  #' Plot TF-gene regulatory network
  #'
  #' @param net Network object from build_pathway_network()
  #' @param title Plot title (defaults to pathway name)
  #' @param layout ggraph layout algorithm

  if (is.null(net)) {
    return(NULL)
  }

  g <- net$graph
  pathway_name <- net$pathway_name

  if (is.null(title)) {
    title <- str_wrap(gsub("_", " ", pathway_name), width = 60)
  }

  # Ensure edge attributes are numeric vectors, not lists
  if (!is.null(E(g)$edge_score)) {
    E(g)$edge_score <- as.numeric(unlist(E(g)$edge_score))
  }
  if (!is.null(E(g)$rna_log2fc)) {
    E(g)$rna_log2fc <- as.numeric(unlist(E(g)$rna_log2fc))
  }

  # Ensure node attributes are numeric
  if (!is.null(V(g)$tf_activity)) {
    V(g)$tf_activity <- as.numeric(unlist(V(g)$tf_activity))
  }
  if (!is.null(V(g)$expression)) {
    V(g)$expression <- as.numeric(unlist(V(g)$expression))
  }

  # Create ggraph
  set.seed(42)  # for reproducibility

  p <- ggraph(g, layout = layout) +
    # Edges - simplified
    geom_edge_link(
      aes(edge_alpha = abs(edge_score)),
      color = "grey50",
      arrow = arrow(length = unit(2, "mm"), type = "closed"),
      end_cap = circle(3, "mm")
    ) +
    scale_edge_alpha_continuous(range = edge_alpha_range, guide = "none") +

    # Gene nodes
    geom_node_point(
      data = . %>% filter(type == "gene"),
      aes(size = degree, fill = expression),
      shape = 21, color = "grey30", stroke = 0.5
    ) +
    scale_size_continuous(range = gene_size_range, name = "Degree") +

    # TF nodes (diamond shape)
    geom_node_point(
      data = . %>% filter(type == "TF"),
      aes(size = degree * 1.5, fill = tf_activity),
      shape = 23, color = "black", stroke = 1
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "Activity/\nExpression",
      na.value = "grey70"
    ) +

    # Labels
    geom_node_text(
      data = . %>% filter(type == "TF"),
      aes(label = name),
      size = 3.5, fontface = "bold",
      repel = TRUE, max.overlaps = 20
    ) +
    geom_node_text(
      data = . %>% filter(type == "gene"),
      aes(label = name),
      size = 2.5,
      repel = TRUE, max.overlaps = 15
    ) +

    # Theme
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      title = title,
      subtitle = sprintf("%d TFs regulating %d genes (%d edges)",
                         net$n_tfs, net$n_genes, net$n_edges)
    )

  return(p)
}

plot_pathway_enrichment_network <- function(edges, nodes, gene_sets,
                                            top_n_pathways = 20,
                                            min_genes = 5) {
  #' Plot emap-style pathway enrichment network
  #'
  #' Shows which pathways share regulated genes

  message("Building pathway enrichment network...")

  # Get all regulated genes
  regulated_genes <- unique(edges$gene)

  # Calculate overlap with each pathway
  pathway_stats <- lapply(names(gene_sets), function(pw) {
    pw_genes <- gene_sets[[pw]]
    overlap <- intersect(regulated_genes, pw_genes)

    data.frame(
      pathway = pw,
      n_overlap = length(overlap),
      n_pathway = length(pw_genes),
      pct_overlap = length(overlap) / length(pw_genes) * 100,
      genes = paste(overlap, collapse = ",")
    )
  }) %>% bind_rows()

  # Filter and rank
  pathway_stats <- pathway_stats %>%
    filter(n_overlap >= min_genes) %>%
    arrange(desc(n_overlap)) %>%
    head(top_n_pathways)

  if (nrow(pathway_stats) < 2) {
    message("  Not enough pathways with sufficient overlap")
    return(NULL)
  }

  # Build pathway-pathway similarity network
  pw_names <- pathway_stats$pathway
  n_pw <- length(pw_names)

  # Calculate Jaccard similarity between pathways (based on regulated genes)
  sim_matrix <- matrix(0, n_pw, n_pw, dimnames = list(pw_names, pw_names))

  for (i in 1:(n_pw-1)) {
    genes_i <- unlist(strsplit(pathway_stats$genes[i], ","))
    for (j in (i+1):n_pw) {
      genes_j <- unlist(strsplit(pathway_stats$genes[j], ","))
      jaccard <- length(intersect(genes_i, genes_j)) /
                 length(union(genes_i, genes_j))
      sim_matrix[i,j] <- jaccard
      sim_matrix[j,i] <- jaccard
    }
  }

  # Create edge list from similarity matrix
  pw_edges <- which(sim_matrix > 0.1, arr.ind = TRUE) %>%
    as.data.frame() %>%
    filter(row < col) %>%
    mutate(
      from = pw_names[row],
      to = pw_names[col],
      similarity = sim_matrix[cbind(row, col)]
    ) %>%
    select(from, to, similarity)

  if (nrow(pw_edges) == 0) {
    message("  No pathway connections found")
    return(NULL)
  }

  # Build graph
  g <- graph_from_data_frame(
    d = pw_edges,
    directed = FALSE,
    vertices = pathway_stats %>% select(pathway, n_overlap, pct_overlap)
  )

  # Clean pathway names for display
  V(g)$label <- sapply(V(g)$name, function(x) {
    x <- gsub("^HALLMARK_|^KEGG_|^REACTOME_|^GO_", "", x)
    x <- gsub("_", " ", x)
    str_wrap(x, width = 20)
  })

  # Plot
  set.seed(42)

  p <- ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(edge_width = similarity, edge_alpha = similarity),
      color = "grey60"
    ) +
    scale_edge_width_continuous(range = c(0.5, 3), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.3, 0.8), guide = "none") +

    geom_node_point(
      aes(size = n_overlap, fill = pct_overlap),
      shape = 21, color = "grey30", stroke = 0.8
    ) +
    scale_size_continuous(range = c(5, 20), name = "Regulated\nGenes") +
    scale_fill_viridis_c(option = "plasma", name = "% Pathway\nRegulated") +

    geom_node_text(
      aes(label = label),
      size = 2.5, repel = TRUE, max.overlaps = 30
    ) +

    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "Pathway Enrichment Network (emap-style)")

  list(
    plot = p,
    pathway_stats = pathway_stats
  )
}

# ==============================================================================
# CYTOSCAPE EXPORT
# ==============================================================================

export_cytoscape <- function(edges, nodes, output_dir, contrast) {
  #' Export network data in Cytoscape-compatible format

  message("Exporting Cytoscape files...")

  # Edge file (SIF-like format for network)
  edges_cy <- edges %>%
    select(tf, gene, edge_score, rna_log2fc, n_peaks) %>%
    mutate(interaction = "regulates")

  write.table(
    edges_cy,
    file.path(output_dir, sprintf("cytoscape_edges_%s.tsv", contrast)),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  # Node attributes file
  nodes_cy <- nodes %>%
    select(node, node_type, everything())

  write.table(
    nodes_cy,
    file.path(output_dir, sprintf("cytoscape_nodes_%s.tsv", contrast)),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  message(sprintf("  Saved to: %s", output_dir))
}

# ==============================================================================
# MAIN
# ==============================================================================

main <- function() {
  # Parse arguments
  parser <- ArgumentParser(description = "TF-Gene Network Visualization")

  parser$add_argument("--edges", type = "character", required = TRUE,
                      help = "Path to edges TSV file")
  parser$add_argument("--nodes", type = "character", required = TRUE,
                      help = "Path to nodes TSV file")
  parser$add_argument("--contrast", type = "character", required = TRUE,
                      help = "Contrast name (for output naming)")
  parser$add_argument("--config", type = "character", default = "config.yaml",
                      help = "Path to config file")
  parser$add_argument("--pathway", type = "character", default = NULL,
                      help = "Specific pathway to plot (e.g., HALLMARK_INTERFERON_ALPHA_RESPONSE)")
  parser$add_argument("--all_pathways", action = "store_true", default = FALSE,
                      help = "Generate plots for all focus pathways")
  parser$add_argument("--emap", action = "store_true", default = FALSE,
                      help = "Generate emap-style pathway enrichment plot")
  parser$add_argument("--output", type = "character", default = NULL,
                      help = "Output directory (overrides config)")

  args <- parser$parse_args()

  # Load config
  cfg <- load_config(args$config)

  # Setup output directory
  if (!is.null(args$output)) {
    out_dir <- args$output
  } else {
    out_dir <- file.path(dirname(args$edges), "..", cfg$output$plots_dir %||% "plots")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Also create cytoscape directory
  cy_dir <- file.path(dirname(out_dir), cfg$output$cytoscape_dir %||% "cytoscape")
  dir.create(cy_dir, recursive = TRUE, showWarnings = FALSE)

  # Load data
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Loading data...")
  message(paste(rep("=", 60), collapse = ""))

  edges <- read.delim(args$edges, stringsAsFactors = FALSE)
  nodes <- read.delim(args$nodes, stringsAsFactors = FALSE)

  message(sprintf("  Loaded %d edges, %d nodes", nrow(edges), nrow(nodes)))

  # Load MSigDB
  gene_sets <- load_msigdb(
    species = cfg$msigdb$species,
    collections = cfg$msigdb$collections
  )

  # Get pathways to plot
  pathways_to_plot <- c()

  if (!is.null(args$pathway)) {
    pathways_to_plot <- args$pathway
  }

  if (args$all_pathways && !is.null(cfg$msigdb$focus_pathways)) {
    pathways_to_plot <- c(pathways_to_plot, cfg$msigdb$focus_pathways)
  }

  # If no specific pathway, just do top pathways based on overlap
  if (length(pathways_to_plot) == 0) {
    # Find pathways with most regulated genes
    regulated_genes <- unique(edges$gene)

    pw_overlap <- sapply(gene_sets, function(g) length(intersect(g, regulated_genes)))
    top_pw <- names(sort(pw_overlap, decreasing = TRUE))[1:min(10, length(pw_overlap))]
    pathways_to_plot <- top_pw[pw_overlap[top_pw] >= 5]

    message(sprintf("\nNo pathway specified, using top %d pathways by gene overlap",
                    length(pathways_to_plot)))
  }

  # Remove duplicates
  pathways_to_plot <- unique(pathways_to_plot)

  # Plot each pathway
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Generating pathway network plots...")
  message(paste(rep("=", 60), collapse = ""))

  for (pw in pathways_to_plot) {
    if (!(pw %in% names(gene_sets))) {
      message(sprintf("  Pathway not found: %s", pw))
      next
    }

    message(sprintf("\nProcessing: %s", pw))

    # Build network
    net <- build_pathway_network(
      edges = edges,
      nodes = nodes,
      pathway_genes = gene_sets[[pw]],
      pathway_name = pw,
      max_genes = cfg$network$max_genes_per_pathway %||% 50,
      max_tfs = cfg$network$max_tfs_per_pathway %||% 15,
      min_tf_activity = cfg$network$min_tf_activity %||% 0.05
    )

    if (is.null(net)) {
      next
    }

    # Create plot
    p <- plot_tf_gene_network(net)

    if (is.null(p)) {
      next
    }

    # Save plot
    base_name <- sprintf("network_%s_%s", args$contrast, safe_name(pw))

    if (cfg$output$save_png %||% TRUE) {
      png_file <- file.path(out_dir, paste0(base_name, ".png"))
      ggsave(png_file, p,
             width = cfg$output$plot_width %||% 12,
             height = cfg$output$plot_height %||% 10,
             dpi = cfg$output$plot_dpi %||% 300)
      message(sprintf("  Saved: %s", basename(png_file)))
    }

    if (cfg$output$save_pdf %||% TRUE) {
      pdf_file <- file.path(out_dir, paste0(base_name, ".pdf"))
      ggsave(pdf_file, p,
             width = cfg$output$plot_width %||% 12,
             height = cfg$output$plot_height %||% 10)
      message(sprintf("  Saved: %s", basename(pdf_file)))
    }
  }

  # Generate emap plot if requested
  if (args$emap) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("Generating emap-style pathway network...")
    message(paste(rep("=", 60), collapse = ""))

    emap_result <- plot_pathway_enrichment_network(
      edges = edges,
      nodes = nodes,
      gene_sets = gene_sets,
      top_n_pathways = 25,
      min_genes = 3
    )

    if (!is.null(emap_result)) {
      base_name <- sprintf("emap_%s", args$contrast)

      if (cfg$output$save_png %||% TRUE) {
        png_file <- file.path(out_dir, paste0(base_name, ".png"))
        ggsave(png_file, emap_result$plot,
               width = 14, height = 12, dpi = 300)
        message(sprintf("  Saved: %s", basename(png_file)))
      }

      if (cfg$output$save_pdf %||% TRUE) {
        pdf_file <- file.path(out_dir, paste0(base_name, ".pdf"))
        ggsave(pdf_file, emap_result$plot,
               width = 14, height = 12)
        message(sprintf("  Saved: %s", basename(pdf_file)))
      }

      # Save pathway stats
      stats_file <- file.path(out_dir, sprintf("pathway_stats_%s.tsv", args$contrast))
      write.table(emap_result$pathway_stats, stats_file,
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }

  # Export Cytoscape files
  if (cfg$output$save_cytoscape %||% TRUE) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("Exporting Cytoscape files...")
    message(paste(rep("=", 60), collapse = ""))

    export_cytoscape(edges, nodes, cy_dir, args$contrast)
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Done!")
  message(paste(rep("=", 60), collapse = ""))
}

# Handle null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Run main
if (!interactive()) {
  main()
}
