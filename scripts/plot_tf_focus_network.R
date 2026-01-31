#!/usr/bin/env Rscript
################################################################################
#  TF-Focused Cascade Network Visualization
#  =========================================
#  Creates hierarchical network plots showing a primary TF with its regulatory
#  cascade - direct targets and secondary TF targets.
#
#  Visual encodings:
#    - Node color: RNA expression log2FC (red=up, blue=down)
#    - Node shape: Diamond=TF, Circle=gene
#    - Node size: Layer-based (L0 largest, L2 smallest)
#    - Edge linetype: Binding strength (dotted→dashed→solid)
#    - Edge width: Differential chromatin accessibility
#
#  Usage:
#    Rscript plot_tf_focus_network.R \
#      --edges tf_focus_edges_STAT1.tsv \
#      --nodes tf_focus_nodes_STAT1.tsv \
#      --tf STAT1 \
#      --contrast "Transformed vs Normal"
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
  library(scales)
  library(stringr)
  library(viridis)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

DEFAULT_CONFIG <- list(
  output = list(
    plot_width = 14,
    plot_height = 12,
    plot_dpi = 300,
    save_pdf = TRUE,
    save_png = TRUE
  ),
  tf_focus = list(
    title_prefix = "TF Regulatory Network",
    contrast_label = "Contrast"
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

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_-]", "_", x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ==============================================================================
# NETWORK BUILDING
# ==============================================================================

build_cascade_graph <- function(edges, nodes) {
  #' Build igraph object from edge/node tables
  #'
  #' @param edges Edge table with source, target, binding_score, accessibility_diff, layer

#' @param nodes Node table with node, node_type, layer, expression

  # Ensure required columns exist
  required_edge_cols <- c("source", "target")
  required_node_cols <- c("node", "node_type", "layer", "expression")

  if (!all(required_edge_cols %in% names(edges))) {
    stop("Missing required edge columns: ", paste(setdiff(required_edge_cols, names(edges)), collapse = ", "))
  }
  if (!all(required_node_cols %in% names(nodes))) {
    stop("Missing required node columns: ", paste(setdiff(required_node_cols, names(nodes)), collapse = ", "))
  }

  # Create graph
  g <- graph_from_data_frame(
    d = edges %>% select(source, target, everything()),
    directed = TRUE,
    vertices = nodes %>% select(node, everything()) %>% rename(name = node)
  )

  return(g)
}

# ==============================================================================
# VISUALIZATION
# ==============================================================================

plot_tf_cascade <- function(g, primary_tf, contrast_label = NULL,
                            layout_type = "sugiyama",
                            output_prefix = NULL,
                            save_pdf = TRUE,
                            save_png = TRUE,
                            width = 14,
                            height = 12,
                            dpi = 300) {
  #' Create cascade network visualization
  #'
  #' @param g igraph object from build_cascade_graph
  #' @param primary_tf Name of primary TF (center of cascade)
  #' @param contrast_label Label for contrast (shown in subtitle)
  #' @param layout_type Layout algorithm: "sugiyama", "tree", "fr", "kk"
  #' @param output_prefix File prefix for saving (without extension)

  message("\nCreating cascade network plot...")

  # Extract node attributes
  node_data <- igraph::as_data_frame(g, what = "vertices")

  # Count nodes per layer
  n_layer0 <- sum(node_data$layer == 0)
  n_layer1 <- sum(node_data$layer == 1)
  n_layer2 <- sum(node_data$layer == 2)
  n_tfs <- sum(node_data$node_type == "TF")
  n_genes <- sum(node_data$node_type == "gene")

  # Secondary TFs (TFs in layer 1)
  secondary_tfs <- node_data %>%
    filter(node_type == "TF", layer == 1) %>%
    pull(name)
  n_secondary <- length(secondary_tfs)

  # Prepare edge attributes
  edge_data <- igraph::as_data_frame(g, what = "edges")

  # Scale binding score for linetype categories
  if ("binding_score_scaled" %in% names(edge_data)) {
    edge_data$binding_cat <- cut(
      edge_data$binding_score_scaled,
      breaks = c(-Inf, 0.33, 0.66, Inf),
      labels = c("weak", "medium", "strong")
    )
  } else {
    edge_data$binding_cat <- "medium"
  }

  # Scale accessibility for edge width
  if ("accessibility_diff_scaled" %in% names(edge_data)) {
    edge_data$acc_width <- rescale(abs(edge_data$accessibility_diff_scaled), to = c(0.3, 2))
  } else {
    edge_data$acc_width <- 1
  }

  # Update graph with processed attributes
  E(g)$binding_cat <- edge_data$binding_cat
  E(g)$acc_width <- edge_data$acc_width

  # Create layout
  set.seed(42)

  # Choose layout based on type
  if (layout_type == "sugiyama") {
    # Sugiyama layout respects hierarchy
    layout <- create_layout(g, layout = "sugiyama")
  } else if (layout_type == "tree") {
    # Tree layout with primary TF as root
    layout <- create_layout(g, layout = "tree", root = which(V(g)$name == primary_tf))
  } else if (layout_type == "fr") {
    # Force-directed layout
    layout <- create_layout(g, layout = "fr", niter = 1000)
  } else if (layout_type == "kk") {
    # Kamada-Kawai
    layout <- create_layout(g, layout = "kk")
  } else {
    # Default to sugiyama
    layout <- create_layout(g, layout = "sugiyama")
  }

  # Define size by layer
  size_by_layer <- c("0" = 20, "1" = 12, "2" = 6)

  # Build plot layer by layer
  p <- ggraph(layout) +

    # --- EDGES ---
    # All edges with binding strength linetype and accessibility width
    geom_edge_link(
      aes(
        edge_linetype = cut(binding_score_scaled,
                            breaks = c(-Inf, 0.33, 0.66, Inf),
                            labels = c("weak", "medium", "strong")),
        edge_width = abs(accessibility_diff_scaled),
        edge_alpha = ifelse(layer == 1, 0.7, 0.4)
      ),
      arrow = arrow(length = unit(2.5, "mm"), type = "closed", angle = 20),
      end_cap = circle(4, "mm"),
      color = "grey40"
    ) +

    # Edge scales
    scale_edge_linetype_manual(
      values = c("weak" = "dotted", "medium" = "dashed", "strong" = "solid"),
      name = "Binding\nStrength",
      guide = guide_legend(order = 2)
    ) +
    scale_edge_width_continuous(
      range = c(0.3, 2.5),
      name = "Chromatin\nOpening",
      guide = guide_legend(order = 3)
    ) +
    scale_edge_alpha_identity() +

    # --- NODES ---
    # Layer 0: Primary TF (largest diamond)
    geom_node_point(
      data = . %>% filter(layer == 0),
      aes(fill = expression),
      shape = 23,  # Diamond
      size = 20,
      stroke = 2.5,
      color = "black"
    ) +

    # Layer 1: Secondary TFs (medium diamonds)
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(fill = expression),
      shape = 23,
      size = 14,
      stroke = 1.8,
      color = "black"
    ) +

    # Layer 1: Terminal genes (medium circles)
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(fill = expression, size = abs(expression) + 4),
      shape = 21,
      stroke = 0.8,
      color = "grey30"
    ) +

    # Layer 2: Downstream genes (small circles)
    geom_node_point(
      data = . %>% filter(layer == 2),
      aes(fill = expression),
      shape = 21,
      size = 5,
      stroke = 0.4,
      color = "grey50"
    ) +

    # --- COLOR SCALE ---
    scale_fill_gradient2(
      low = "#2166AC",   # Blue (down-regulated)
      mid = "white",
      high = "#B2182B",  # Red (up-regulated)
      midpoint = 0,
      name = "Expression\nlog2FC",
      guide = guide_colorbar(order = 1)
    ) +

    # Size scale (only for variable-sized nodes)
    scale_size_continuous(
      range = c(6, 12),
      guide = "none"  # Hide size legend (explained in caption)
    ) +

    # --- LABELS ---
    # Primary TF label (bold, large, centered on node)
    geom_node_text(
      data = . %>% filter(layer == 0),
      aes(label = name),
      size = 5,
      fontface = "bold",
      vjust = 0.5,
      color = "black"
    ) +

    # Secondary TF labels (bold, medium, centered on node)
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(label = name),
      size = 3.5,
      fontface = "bold",
      vjust = 0.5,
      color = "black"
    ) +

    # Layer 1 gene labels (normal, centered on node)
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(label = name),
      size = 2.5,
      vjust = 0.5,
      color = "black"
    ) +

    # Layer 2 labels (small, centered on node)
    geom_node_text(
      data = . %>% filter(layer == 2),
      aes(label = name),
      size = 2,
      color = "black",
      vjust = 0.5
    ) +

    # --- THEME ---
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30",
                                   margin = margin(b = 10)),
      plot.caption = element_text(hjust = 0, size = 9, color = "grey50",
                                  margin = margin(t = 15)),
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin(l = 10),
      plot.margin = margin(15, 15, 15, 15)
    ) +

    # --- TITLES ---
    labs(
      title = paste(primary_tf, "Regulatory Cascade"),
      subtitle = paste0(
        "Contrast: ", contrast_label %||% "Contrast", "\n",
        "Layer 1: ", n_layer1, " direct targets | ",
        "Secondary TFs: ", n_secondary,
        ifelse(n_secondary > 0, paste0(" (", paste(secondary_tfs, collapse = ", "), ")"), ""),
        "\n",
        "Layer 2: ", n_layer2, " downstream targets"
      ),
      caption = paste0(
        "Node shapes: ", intToUtf8(9670), " = TF, ", intToUtf8(9679), " = Gene\n",
        "Node size indicates regulatory layer (larger = closer to ", primary_tf, ")\n",
        "Edge style: dotted = weak binding, solid = strong binding\n",
        "Edge width: thicker = greater chromatin accessibility change"
      )
    )

  # Save plot
  if (!is.null(output_prefix)) {
    if (save_png) {
      png_file <- paste0(output_prefix, ".png")
      ggsave(png_file, p, width = width, height = height, dpi = dpi)
      message(sprintf("  Saved: %s", basename(png_file)))
    }

    if (save_pdf) {
      pdf_file <- paste0(output_prefix, ".pdf")
      ggsave(pdf_file, p, width = width, height = height)
      message(sprintf("  Saved: %s", basename(pdf_file)))
    }
  }

  return(p)
}

plot_tf_cascade_circular <- function(g, primary_tf, contrast_label = NULL,
                                     output_prefix = NULL,
                                     save_pdf = TRUE,
                                     save_png = TRUE,
                                     width = 14,
                                     height = 14,
                                     dpi = 300) {
  #' Alternative circular/radial layout for cascade network
  #'
  #' Primary TF at center, Layer 1 in inner ring, Layer 2 in outer ring

  message("\nCreating circular cascade network plot...")

  # Extract node attributes
  node_data <- igraph::as_data_frame(g, what = "vertices")
  edge_data <- igraph::as_data_frame(g, what = "edges")

  # Count nodes
  n_layer1 <- sum(node_data$layer == 1)
  n_layer2 <- sum(node_data$layer == 2)
  secondary_tfs <- node_data %>%
    filter(node_type == "TF", layer == 1) %>%
    pull(name)
  n_secondary <- length(secondary_tfs)

  # Scale binding score for linetype
  if ("binding_score_scaled" %in% names(edge_data)) {
    edge_data$binding_cat <- cut(
      edge_data$binding_score_scaled,
      breaks = c(-Inf, 0.33, 0.66, Inf),
      labels = c("weak", "medium", "strong")
    )
  } else {
    edge_data$binding_cat <- "medium"
  }

  # Scale accessibility
  if ("accessibility_diff_scaled" %in% names(edge_data)) {
    edge_data$acc_width <- rescale(abs(edge_data$accessibility_diff_scaled), to = c(0.3, 2))
  } else {
    edge_data$acc_width <- 1
  }

  # Update graph
  E(g)$binding_cat <- edge_data$binding_cat
  E(g)$acc_width <- edge_data$acc_width

  # Create manual radial layout
  # Primary TF at center (0,0)
  # Layer 1 in inner ring (radius 1)
  # Layer 2 in outer ring (radius 2)

  node_data <- node_data %>%
    group_by(layer) %>%
    mutate(
      idx = row_number(),
      n_in_layer = n(),
      angle = 2 * pi * (idx - 1) / n_in_layer
    ) %>%
    ungroup() %>%
    mutate(
      radius = case_when(
        layer == 0 ~ 0,
        layer == 1 ~ 1,
        layer == 2 ~ 2.2
      ),
      x = radius * cos(angle),
      y = radius * sin(angle)
    )

  # Create layout data frame
  layout_df <- node_data %>%
    select(name, x, y) %>%
    as.data.frame()

  # Match order to graph vertices
  vertex_order <- V(g)$name
  layout_matrix <- layout_df %>%
    slice(match(vertex_order, name)) %>%
    select(x, y) %>%
    as.matrix()

  # Create ggraph layout
  layout <- create_layout(g, layout = "manual", x = layout_matrix[,1], y = layout_matrix[,2])

  # Build plot
  p <- ggraph(layout) +

    # Add concentric circles as guides
    annotate("path",
             x = cos(seq(0, 2*pi, length.out = 100)),
             y = sin(seq(0, 2*pi, length.out = 100)),
             color = "grey90", linetype = "dashed", linewidth = 0.3) +
    annotate("path",
             x = 2.2 * cos(seq(0, 2*pi, length.out = 100)),
             y = 2.2 * sin(seq(0, 2*pi, length.out = 100)),
             color = "grey90", linetype = "dashed", linewidth = 0.3) +

    # Edges
    geom_edge_link(
      aes(
        edge_linetype = cut(binding_score_scaled,
                            breaks = c(-Inf, 0.33, 0.66, Inf),
                            labels = c("weak", "medium", "strong")),
        edge_width = abs(accessibility_diff_scaled),
        edge_alpha = ifelse(layer == 1, 0.7, 0.4)
      ),
      arrow = arrow(length = unit(2, "mm"), type = "closed", angle = 20),
      end_cap = circle(3, "mm"),
      color = "grey40"
    ) +

    scale_edge_linetype_manual(
      values = c("weak" = "dotted", "medium" = "dashed", "strong" = "solid"),
      name = "Binding\nStrength"
    ) +
    scale_edge_width_continuous(range = c(0.3, 2), name = "Chromatin\nOpening") +
    scale_edge_alpha_identity() +

    # Nodes
    geom_node_point(
      data = . %>% filter(layer == 0),
      aes(fill = expression),
      shape = 23, size = 22, stroke = 3, color = "black"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(fill = expression),
      shape = 23, size = 14, stroke = 2, color = "black"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(fill = expression),
      shape = 21, size = 10, stroke = 1, color = "grey30"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 2),
      aes(fill = expression),
      shape = 21, size = 5, stroke = 0.5, color = "grey50"
    ) +

    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "Expression\nlog2FC"
    ) +

    # Labels - centered on nodes
    geom_node_text(
      data = . %>% filter(layer == 0),
      aes(label = name),
      size = 5, fontface = "bold", vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(label = name),
      size = 3.5, fontface = "bold", vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(label = name),
      size = 2.5, vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 2),
      aes(label = name),
      size = 2, color = "black", vjust = 0.5
    ) +

    coord_fixed() +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
      plot.caption = element_text(hjust = 0, size = 9, color = "grey50"),
      legend.position = "right",
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(
      title = paste(primary_tf, "Regulatory Cascade"),
      subtitle = paste0(
        "Contrast: ", contrast_label %||% "Contrast", "\n",
        "Inner ring: ", n_layer1, " direct targets | ",
        "Secondary TFs: ", n_secondary, "\n",
        "Outer ring: ", n_layer2, " downstream targets"
      ),
      caption = "Node shapes: diamond = TF, circle = gene | Size indicates layer"
    )

  # Save
  if (!is.null(output_prefix)) {
    if (save_png) {
      png_file <- paste0(output_prefix, "_circular.png")
      ggsave(png_file, p, width = width, height = height, dpi = dpi)
      message(sprintf("  Saved: %s", basename(png_file)))
    }
    if (save_pdf) {
      pdf_file <- paste0(output_prefix, "_circular.pdf")
      ggsave(pdf_file, p, width = width, height = height)
      message(sprintf("  Saved: %s", basename(pdf_file)))
    }
  }

  return(p)
}

# ==============================================================================
# COMPACT PLOTS (for figure panels)
# ==============================================================================

plot_tf_cascade_compact <- function(g, primary_tf, contrast_label = NULL,
                                    output_prefix = NULL,
                                    save_pdf = TRUE,
                                    save_png = TRUE,
                                    width = 4,
                                    height = 3.5,
                                    dpi = 300) {
  #' Create compact cascade network for figure panels
  #'
  #' Smaller size (4" x 3.5"), minimal legends, smaller text

  message("\nCreating compact cascade network plot...")

  # Extract node attributes
  node_data <- igraph::as_data_frame(g, what = "vertices")
  edge_data <- igraph::as_data_frame(g, what = "edges")

  # Count nodes
  n_layer1 <- sum(node_data$layer == 1)
  n_layer2 <- sum(node_data$layer == 2)
  secondary_tfs <- node_data %>%
    filter(node_type == "TF", layer == 1) %>%
    pull(name)
  n_secondary <- length(secondary_tfs)

  # Scale binding score for linetype
  if ("binding_score_scaled" %in% names(edge_data)) {
    edge_data$binding_cat <- cut(
      edge_data$binding_score_scaled,
      breaks = c(-Inf, 0.33, 0.66, Inf),
      labels = c("weak", "medium", "strong")
    )
  } else {
    edge_data$binding_cat <- "medium"
  }

  # Update graph
  E(g)$binding_cat <- edge_data$binding_cat

  # Create layout
  set.seed(42)
  layout <- create_layout(g, layout = "sugiyama")

  # Build compact plot
  p <- ggraph(layout) +

    # Edges - simplified
    geom_edge_link(
      aes(
        edge_linetype = cut(binding_score_scaled,
                            breaks = c(-Inf, 0.33, 0.66, Inf),
                            labels = c("weak", "medium", "strong")),
        edge_alpha = ifelse(layer == 1, 0.7, 0.4)
      ),
      edge_width = 0.4,
      arrow = arrow(length = unit(1.5, "mm"), type = "closed", angle = 20),
      end_cap = circle(2, "mm"),
      color = "grey40"
    ) +

    scale_edge_linetype_manual(
      values = c("weak" = "dotted", "medium" = "dashed", "strong" = "solid"),
      guide = "none"
    ) +
    scale_edge_alpha_identity() +

    # Nodes - smaller sizes
    geom_node_point(
      data = . %>% filter(layer == 0),
      aes(fill = expression),
      shape = 23, size = 10, stroke = 1.5, color = "black"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(fill = expression),
      shape = 23, size = 7, stroke = 1, color = "black"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(fill = expression),
      shape = 21, size = 5, stroke = 0.5, color = "grey30"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 2),
      aes(fill = expression),
      shape = 21, size = 3, stroke = 0.3, color = "grey50"
    ) +

    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "log2FC",
      guide = guide_colorbar(barwidth = 0.5, barheight = 3)
    ) +

    # Labels - smaller text
    geom_node_text(
      data = . %>% filter(layer == 0),
      aes(label = name),
      size = 2.5, fontface = "bold", vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(label = name),
      size = 2, fontface = "bold", vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(label = name),
      size = 1.5, vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 2),
      aes(label = name),
      size = 1.2, color = "black", vjust = 0.5
    ) +

    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 7, color = "grey30"),
      legend.position = "right",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    labs(
      title = paste(primary_tf, "Network"),
      subtitle = paste0("L1: ", n_layer1, " | L2: ", n_layer2)
    )

  # Save
  if (!is.null(output_prefix)) {
    if (save_png) {
      png_file <- paste0(output_prefix, "_compact.png")
      ggsave(png_file, p, width = width, height = height, dpi = dpi)
      message(sprintf("  Saved: %s", basename(png_file)))
    }
    if (save_pdf) {
      pdf_file <- paste0(output_prefix, "_compact.pdf")
      ggsave(pdf_file, p, width = width, height = height)
      message(sprintf("  Saved: %s", basename(pdf_file)))
    }
  }

  return(p)
}

plot_tf_cascade_circular_compact <- function(g, primary_tf, contrast_label = NULL,
                                              output_prefix = NULL,
                                              save_pdf = TRUE,
                                              save_png = TRUE,
                                              width = 4,
                                              height = 4,
                                              dpi = 300) {
  #' Create compact circular cascade network for figure panels
  #'
  #' Smaller size (4" x 4"), minimal legends, smaller text

  message("\nCreating compact circular cascade network plot...")

  # Extract node attributes
  node_data <- igraph::as_data_frame(g, what = "vertices")
  edge_data <- igraph::as_data_frame(g, what = "edges")

  # Count nodes
  n_layer1 <- sum(node_data$layer == 1)
  n_layer2 <- sum(node_data$layer == 2)

  # Scale binding score for linetype
  if ("binding_score_scaled" %in% names(edge_data)) {
    edge_data$binding_cat <- cut(
      edge_data$binding_score_scaled,
      breaks = c(-Inf, 0.33, 0.66, Inf),
      labels = c("weak", "medium", "strong")
    )
  } else {
    edge_data$binding_cat <- "medium"
  }

  # Update graph
  E(g)$binding_cat <- edge_data$binding_cat

  # Create manual radial layout
  node_data <- node_data %>%
    group_by(layer) %>%
    mutate(
      idx = row_number(),
      n_in_layer = n(),
      angle = 2 * pi * (idx - 1) / n_in_layer
    ) %>%
    ungroup() %>%
    mutate(
      radius = case_when(
        layer == 0 ~ 0,
        layer == 1 ~ 1,
        layer == 2 ~ 2.2
      ),
      x = radius * cos(angle),
      y = radius * sin(angle)
    )

  # Create layout data frame
  layout_df <- node_data %>%
    select(name, x, y) %>%
    as.data.frame()

  # Match order to graph vertices
  vertex_order <- V(g)$name
  layout_matrix <- layout_df %>%
    slice(match(vertex_order, name)) %>%
    select(x, y) %>%
    as.matrix()

  # Create ggraph layout
  layout <- create_layout(g, layout = "manual", x = layout_matrix[,1], y = layout_matrix[,2])

  # Build compact circular plot
  p <- ggraph(layout) +

    # Concentric circle guides
    annotate("path",
             x = cos(seq(0, 2*pi, length.out = 100)),
             y = sin(seq(0, 2*pi, length.out = 100)),
             color = "grey90", linetype = "dashed", linewidth = 0.2) +
    annotate("path",
             x = 2.2 * cos(seq(0, 2*pi, length.out = 100)),
             y = 2.2 * sin(seq(0, 2*pi, length.out = 100)),
             color = "grey90", linetype = "dashed", linewidth = 0.2) +

    # Edges - simplified
    geom_edge_link(
      aes(
        edge_linetype = cut(binding_score_scaled,
                            breaks = c(-Inf, 0.33, 0.66, Inf),
                            labels = c("weak", "medium", "strong")),
        edge_alpha = ifelse(layer == 1, 0.6, 0.3)
      ),
      edge_width = 0.3,
      arrow = arrow(length = unit(1, "mm"), type = "closed", angle = 20),
      end_cap = circle(1.5, "mm"),
      color = "grey40"
    ) +

    scale_edge_linetype_manual(
      values = c("weak" = "dotted", "medium" = "dashed", "strong" = "solid"),
      guide = "none"
    ) +
    scale_edge_alpha_identity() +

    # Nodes - smaller sizes
    geom_node_point(
      data = . %>% filter(layer == 0),
      aes(fill = expression),
      shape = 23, size = 8, stroke = 1.2, color = "black"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(fill = expression),
      shape = 23, size = 5, stroke = 0.8, color = "black"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(fill = expression),
      shape = 21, size = 4, stroke = 0.4, color = "grey30"
    ) +
    geom_node_point(
      data = . %>% filter(layer == 2),
      aes(fill = expression),
      shape = 21, size = 2.5, stroke = 0.2, color = "grey50"
    ) +

    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "log2FC",
      guide = guide_colorbar(barwidth = 0.4, barheight = 2.5)
    ) +

    # Labels - smaller text
    geom_node_text(
      data = . %>% filter(layer == 0),
      aes(label = name),
      size = 2, fontface = "bold", vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(label = name),
      size = 1.5, fontface = "bold", vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(label = name),
      size = 1.2, vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 2),
      aes(label = name),
      size = 1, color = "black", vjust = 0.5
    ) +

    coord_fixed() +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 7, color = "grey30"),
      legend.position = "right",
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 4),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    labs(
      title = paste(primary_tf, "Network"),
      subtitle = paste0("L1: ", n_layer1, " | L2: ", n_layer2)
    )

  # Save
  if (!is.null(output_prefix)) {
    if (save_png) {
      png_file <- paste0(output_prefix, "_circular_compact.png")
      ggsave(png_file, p, width = width, height = height, dpi = dpi)
      message(sprintf("  Saved: %s", basename(png_file)))
    }
    if (save_pdf) {
      pdf_file <- paste0(output_prefix, "_circular_compact.pdf")
      ggsave(pdf_file, p, width = width, height = height)
      message(sprintf("  Saved: %s", basename(pdf_file)))
    }
  }

  return(p)
}

# ==============================================================================
# MAIN
# ==============================================================================

main <- function() {
  # Parse arguments
  parser <- ArgumentParser(description = "TF-Focused Cascade Network Visualization")

  parser$add_argument("--edges", type = "character", required = TRUE,
                      help = "Path to edges TSV file from build_tf_focus_network.py")
  parser$add_argument("--nodes", type = "character", required = TRUE,
                      help = "Path to nodes TSV file")
  parser$add_argument("--tf", type = "character", required = TRUE,
                      help = "Primary TF name (e.g., STAT1)")
  parser$add_argument("--contrast", type = "character", default = NULL,
                      help = "Contrast label for subtitle")
  parser$add_argument("--config", type = "character", default = "config.yaml",
                      help = "Path to config file")
  parser$add_argument("--output", type = "character", default = NULL,
                      help = "Output directory (default: same as edges file)")
  parser$add_argument("--layout", type = "character", default = "sugiyama",
                      help = "Layout type: sugiyama, tree, fr, kk, circular")
  parser$add_argument("--circular", action = "store_true", default = FALSE,
                      help = "Also generate circular layout version")
  parser$add_argument("--compact", action = "store_true", default = FALSE,
                      help = "Also generate compact version for figure panels (4x3.5 inches)")

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

  # Load data
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Loading data...")
  message(paste(rep("=", 60), collapse = ""))

  edges <- read.delim(args$edges, stringsAsFactors = FALSE)
  nodes <- read.delim(args$nodes, stringsAsFactors = FALSE)

  message(sprintf("  Loaded %d edges, %d nodes", nrow(edges), nrow(nodes)))

  # Validate data
  if (nrow(edges) == 0 || nrow(nodes) == 0) {
    stop("No data in edges or nodes files!")
  }

  # Build graph
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Building network graph...")
  message(paste(rep("=", 60), collapse = ""))

  g <- build_cascade_graph(edges, nodes)

  message(sprintf("  Graph: %d vertices, %d edges", vcount(g), ecount(g)))

  # Get contrast label
  contrast_label <- args$contrast
  if (is.null(contrast_label) && !is.null(cfg$tf_focus$contrast_label)) {
    contrast_label <- cfg$tf_focus$contrast_label
  }

  # Create plots
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Generating plots...")
  message(paste(rep("=", 60), collapse = ""))

  output_prefix <- file.path(out_dir, sprintf("tf_focus_%s", safe_name(args$tf)))

  # Main plot
  p1 <- plot_tf_cascade(
    g = g,
    primary_tf = args$tf,
    contrast_label = contrast_label,
    layout_type = args$layout,
    output_prefix = output_prefix,
    save_pdf = cfg$output$save_pdf %||% TRUE,
    save_png = cfg$output$save_png %||% TRUE,
    width = cfg$output$plot_width %||% 14,
    height = cfg$output$plot_height %||% 12,
    dpi = cfg$output$plot_dpi %||% 300
  )

  # Circular version if requested
  if (args$circular || args$layout == "circular") {
    p2 <- plot_tf_cascade_circular(
      g = g,
      primary_tf = args$tf,
      contrast_label = contrast_label,
      output_prefix = output_prefix,
      save_pdf = cfg$output$save_pdf %||% TRUE,
      save_png = cfg$output$save_png %||% TRUE,
      width = 14,
      height = 14,
      dpi = cfg$output$plot_dpi %||% 300
    )
  }

  # Compact version if requested
  if (args$compact) {
    p3 <- plot_tf_cascade_compact(
      g = g,
      primary_tf = args$tf,
      contrast_label = contrast_label,
      output_prefix = output_prefix,
      save_pdf = cfg$output$save_pdf %||% TRUE,
      save_png = cfg$output$save_png %||% TRUE,
      width = 4,
      height = 3.5,
      dpi = cfg$output$plot_dpi %||% 300
    )

    # Compact circular version if both --compact and --circular are set
    if (args$circular || args$layout == "circular") {
      p4 <- plot_tf_cascade_circular_compact(
        g = g,
        primary_tf = args$tf,
        contrast_label = contrast_label,
        output_prefix = output_prefix,
        save_pdf = cfg$output$save_pdf %||% TRUE,
        save_png = cfg$output$save_png %||% TRUE,
        width = 4,
        height = 4,
        dpi = cfg$output$plot_dpi %||% 300
      )
    }
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Done!")
  message(paste(rep("=", 60), collapse = ""))
}

# Run main
if (!interactive()) {
  main()
}
