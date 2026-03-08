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
#    - Edge linetype: Binding strength (dotted->dashed->solid)
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
# SHARED VISUALIZATION HELPERS
# ==============================================================================

prepare_edge_attributes <- function(g) {
  #' Compute binding category and accessibility width from graph edge data.
  #' Updates graph edge attributes in place and returns the modified graph.
  #'
  #' @param g igraph object
  #' @return igraph object with binding_cat and acc_width edge attributes

  edge_data <- igraph::as_data_frame(g, what = "edges")

  # Scale binding score into linetype categories

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

  # Update graph
  E(g)$binding_cat <- edge_data$binding_cat
  E(g)$acc_width <- edge_data$acc_width

  return(g)
}

build_node_summary <- function(g) {
  #' Extract node layer counts and secondary TF info from graph.
  #'
  #' @param g igraph object
  #' @return list with n_layer0, n_layer1, n_layer2, n_tfs, n_genes,
  #'         secondary_tfs, n_secondary

  node_data <- igraph::as_data_frame(g, what = "vertices")

  secondary_tfs <- node_data %>%
    filter(node_type == "TF", layer == 1) %>%
    pull(name)

  list(
    n_layer0    = sum(node_data$layer == 0),
    n_layer1    = sum(node_data$layer == 1),
    n_layer2    = sum(node_data$layer == 2),
    n_tfs       = sum(node_data$node_type == "TF"),
    n_genes     = sum(node_data$node_type == "gene"),
    secondary_tfs = secondary_tfs,
    n_secondary   = length(secondary_tfs)
  )
}

build_radial_layout <- function(g) {
  #' Create a manual radial layout for the graph.
  #' Primary TF at center (0,0), Layer 1 at radius 1, Layer 2 at radius 2.2.
  #'
  #' @param g igraph object
  #' @return ggraph layout object

  node_data <- igraph::as_data_frame(g, what = "vertices")

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

  layout_df <- node_data %>%
    select(name, x, y) %>%
    as.data.frame()

  vertex_order <- V(g)$name
  layout_matrix <- layout_df %>%
    slice(match(vertex_order, name)) %>%
    select(x, y) %>%
    as.matrix()

  create_layout(g, layout = "manual", x = layout_matrix[,1], y = layout_matrix[,2])
}

apply_cascade_theme <- function(p, compact = FALSE) {
  #' Apply shared ggplot theme elements for cascade plots.
  #'
  #' @param p ggplot object
  #' @param compact logical; if TRUE use smaller text and margins
  #' @return ggplot object with theme applied

  if (compact) {
    p <- p +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 7, color = "grey30"),
        legend.position = "right",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        plot.margin = margin(5, 5, 5, 5)
      )
  } else {
    p <- p +
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
      )
  }

  return(p)
}

save_cascade_plot <- function(p, output_prefix, suffix = "",
                              save_pdf = TRUE, save_png = TRUE,
                              width = 14, height = 12, dpi = 300) {
  #' Save a cascade plot to PNG and/or PDF.
  #'
  #' @param p ggplot object
  #' @param output_prefix file path prefix (no extension)
  #' @param suffix string appended before extension (e.g. "_circular", "_compact")
  #' @param save_pdf logical
  #' @param save_png logical
  #' @param width plot width in inches
  #' @param height plot height in inches
  #' @param dpi resolution for PNG

  if (!is.null(output_prefix)) {
    if (save_png) {
      png_file <- paste0(output_prefix, suffix, ".png")
      ggsave(png_file, p, width = width, height = height, dpi = dpi)
      message(sprintf("  Saved: %s", basename(png_file)))
    }
    if (save_pdf) {
      pdf_file <- paste0(output_prefix, suffix, ".pdf")
      ggsave(pdf_file, p, width = width, height = height)
      message(sprintf("  Saved: %s", basename(pdf_file)))
    }
  }
}

# ==============================================================================
# CORE VISUALIZATION FUNCTION
# ==============================================================================

plot_tf_cascade_core <- function(g, primary_tf, contrast_label = NULL,
                                  layout_type = "sugiyama",
                                  output_prefix = NULL,
                                  save_pdf = TRUE,
                                  save_png = TRUE,
                                  width = 14,
                                  height = 12,
                                  dpi = 300,
                                  circular = FALSE,
                                  compact = FALSE) {
  #' Core cascade network visualization handling all 4 variants.
  #'
  #' @param g igraph object from build_cascade_graph
  #' @param primary_tf Name of primary TF (center of cascade)
  #' @param contrast_label Label for contrast (shown in subtitle)
  #' @param layout_type Layout algorithm: "sugiyama", "tree", "fr", "kk"
  #'        (ignored when circular = TRUE)
  #' @param output_prefix File prefix for saving (without extension)
  #' @param save_pdf logical
  #' @param save_png logical
  #' @param width plot width in inches
  #' @param height plot height in inches
  #' @param dpi resolution for PNG
  #' @param circular logical; if TRUE use radial layout
  #' @param compact logical; if TRUE use smaller sizes for figure panels
  #' @return ggplot object

  variant_label <- paste0(
    ifelse(compact, "compact ", ""),
    ifelse(circular, "circular ", ""),
    "cascade network plot"
  )
  message(sprintf("\nCreating %s...", variant_label))

  # --- Prepare graph edge attributes ---
  g <- prepare_edge_attributes(g)

  # --- Node summary stats ---
  ns <- build_node_summary(g)

  # --- Sizing parameters (compact vs full) ---
  if (compact) {
    # Node sizes
    sz_l0 <- 10; stroke_l0 <- 1.5
    sz_l1_tf <- 7; stroke_l1_tf <- 1
    sz_l1_gene <- 5; stroke_l1_gene <- 0.5
    sz_l2 <- 3; stroke_l2 <- 0.3
    # Label sizes
    lbl_l0 <- 2.5; lbl_l1_tf <- 2; lbl_l1_gene <- 1.5; lbl_l2 <- 1.2
    # Edge parameters
    arrow_len <- unit(1.5, "mm"); end_cap_r <- 2
    edge_alpha_l1 <- 0.7; edge_alpha_other <- 0.4
    fixed_edge_width <- 0.4  # compact uses fixed width
    # Colorbar
    colorbar_guide <- guide_colorbar(barwidth = 0.5, barheight = 3)
    colorbar_name <- "log2FC"
  } else {
    sz_l0 <- 20; stroke_l0 <- 2.5
    sz_l1_tf <- 14; stroke_l1_tf <- 1.8
    sz_l1_gene <- 10; stroke_l1_gene <- 1
    sz_l2 <- 5; stroke_l2 <- 0.5
    lbl_l0 <- 5; lbl_l1_tf <- 3.5; lbl_l1_gene <- 2.5; lbl_l2 <- 2
    arrow_len <- unit(2.5, "mm"); end_cap_r <- 4
    edge_alpha_l1 <- 0.7; edge_alpha_other <- 0.4
    fixed_edge_width <- NULL  # full uses accessibility-scaled width
    colorbar_guide <- guide_colorbar(order = 1)
    colorbar_name <- "Expression\nlog2FC"
  }

  # Circular variants use slightly different sizes
  if (circular && !compact) {
    sz_l0 <- 22; stroke_l0 <- 3
    sz_l1_tf <- 14; stroke_l1_tf <- 2
    sz_l1_gene <- 10; stroke_l1_gene <- 1
    sz_l2 <- 5; stroke_l2 <- 0.5
    arrow_len <- unit(2, "mm"); end_cap_r <- 3
    colorbar_guide <- guide_colorbar()
    colorbar_name <- "Expression\nlog2FC"
  } else if (circular && compact) {
    sz_l0 <- 8; stroke_l0 <- 1.2
    sz_l1_tf <- 5; stroke_l1_tf <- 0.8
    sz_l1_gene <- 4; stroke_l1_gene <- 0.4
    sz_l2 <- 2.5; stroke_l2 <- 0.2
    lbl_l0 <- 2; lbl_l1_tf <- 1.5; lbl_l1_gene <- 1.2; lbl_l2 <- 1
    arrow_len <- unit(1, "mm"); end_cap_r <- 1.5
    edge_alpha_l1 <- 0.6; edge_alpha_other <- 0.3
    fixed_edge_width <- 0.3
    colorbar_guide <- guide_colorbar(barwidth = 0.4, barheight = 2.5)
  }

  # --- Layout ---
  set.seed(42)

  if (circular) {
    layout <- build_radial_layout(g)
  } else {
    if (layout_type == "tree") {
      layout <- create_layout(g, layout = "tree", root = which(V(g)$name == primary_tf))
    } else if (layout_type == "fr") {
      layout <- create_layout(g, layout = "fr", niter = 1000)
    } else if (layout_type == "kk") {
      layout <- create_layout(g, layout = "kk")
    } else {
      # Default: sugiyama
      layout <- create_layout(g, layout = "sugiyama")
    }
  }

  # --- Build plot ---
  p <- ggraph(layout)

  # Concentric circle guides for circular layouts
  if (circular) {
    guide_lw <- ifelse(compact, 0.2, 0.3)
    p <- p +
      annotate("path",
               x = cos(seq(0, 2*pi, length.out = 100)),
               y = sin(seq(0, 2*pi, length.out = 100)),
               color = "grey90", linetype = "dashed", linewidth = guide_lw) +
      annotate("path",
               x = 2.2 * cos(seq(0, 2*pi, length.out = 100)),
               y = 2.2 * sin(seq(0, 2*pi, length.out = 100)),
               color = "grey90", linetype = "dashed", linewidth = guide_lw)
  }

  # --- Edges ---
  if (!is.null(fixed_edge_width)) {
    # Compact: fixed edge width, no accessibility scaling
    p <- p +
      geom_edge_link(
        aes(
          edge_linetype = cut(binding_score_scaled,
                              breaks = c(-Inf, 0.33, 0.66, Inf),
                              labels = c("weak", "medium", "strong")),
          edge_alpha = ifelse(layer == 1, edge_alpha_l1, edge_alpha_other)
        ),
        edge_width = fixed_edge_width,
        arrow = arrow(length = arrow_len, type = "closed", angle = 20),
        end_cap = circle(end_cap_r, "mm"),
        color = "grey40"
      )
  } else {
    # Full: accessibility-scaled edge width
    p <- p +
      geom_edge_link(
        aes(
          edge_linetype = cut(binding_score_scaled,
                              breaks = c(-Inf, 0.33, 0.66, Inf),
                              labels = c("weak", "medium", "strong")),
          edge_width = abs(accessibility_diff_scaled),
          edge_alpha = ifelse(layer == 1, edge_alpha_l1, edge_alpha_other)
        ),
        arrow = arrow(length = arrow_len, type = "closed", angle = 20),
        end_cap = circle(end_cap_r, "mm"),
        color = "grey40"
      )
  }

  # Edge scales
  if (compact) {
    p <- p +
      scale_edge_linetype_manual(
        values = c("weak" = "dotted", "medium" = "dashed", "strong" = "solid"),
        guide = "none"
      ) +
      scale_edge_alpha_identity()
  } else {
    edge_lt_guide <- guide_legend(order = 2)
    p <- p +
      scale_edge_linetype_manual(
        values = c("weak" = "dotted", "medium" = "dashed", "strong" = "solid"),
        name = "Binding\nStrength",
        guide = edge_lt_guide
      ) +
      scale_edge_alpha_identity()

    if (is.null(fixed_edge_width)) {
      p <- p +
        scale_edge_width_continuous(
          range = c(0.3, 2.5),
          name = "Chromatin\nOpening",
          guide = guide_legend(order = 3)
        )
    }
  }

  # --- Nodes ---
  # Layer 0: Primary TF (diamond)
  p <- p +
    geom_node_point(
      data = . %>% filter(layer == 0),
      aes(fill = expression),
      shape = 23, size = sz_l0, stroke = stroke_l0, color = "black"
    )

  # Layer 1: Secondary TFs (diamond)
  p <- p +
    geom_node_point(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(fill = expression),
      shape = 23, size = sz_l1_tf, stroke = stroke_l1_tf, color = "black"
    )

  # Layer 1: Terminal genes (circle)
  # Full non-circular uses expression-based sizing; all others use fixed size
  if (!compact && !circular) {
    p <- p +
      geom_node_point(
        data = . %>% filter(layer == 1, node_type == "gene"),
        aes(fill = expression, size = abs(expression) + 4),
        shape = 21, stroke = 0.8, color = "grey30"
      ) +
      scale_size_continuous(
        range = c(6, 12),
        guide = "none"
      )
  } else {
    p <- p +
      geom_node_point(
        data = . %>% filter(layer == 1, node_type == "gene"),
        aes(fill = expression),
        shape = 21, size = sz_l1_gene, stroke = stroke_l1_gene, color = "grey30"
      )
  }

  # Layer 2: Downstream genes (circle)
  p <- p +
    geom_node_point(
      data = . %>% filter(layer == 2),
      aes(fill = expression),
      shape = 21, size = sz_l2, stroke = stroke_l2, color = "grey50"
    )

  # --- Color scale ---
  p <- p +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = colorbar_name,
      guide = colorbar_guide
    )

  # --- Labels ---
  p <- p +
    geom_node_text(
      data = . %>% filter(layer == 0),
      aes(label = name),
      size = lbl_l0, fontface = "bold", vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "TF"),
      aes(label = name),
      size = lbl_l1_tf, fontface = "bold", vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 1, node_type == "gene"),
      aes(label = name),
      size = lbl_l1_gene, vjust = 0.5, color = "black"
    ) +
    geom_node_text(
      data = . %>% filter(layer == 2),
      aes(label = name),
      size = lbl_l2, color = "black", vjust = 0.5
    )

  # --- coord_fixed for circular layouts ---
  if (circular) {
    p <- p + coord_fixed()
  }

  # --- Theme ---
  p <- apply_cascade_theme(p, compact = compact)

  # --- Titles ---
  if (compact) {
    p <- p +
      labs(
        title = paste(primary_tf, "Network"),
        subtitle = paste0("L1: ", ns$n_layer1, " | L2: ", ns$n_layer2)
      )
  } else if (circular) {
    p <- p +
      labs(
        title = paste(primary_tf, "Regulatory Cascade"),
        subtitle = paste0(
          "Contrast: ", contrast_label %||% "Contrast", "\n",
          "Inner ring: ", ns$n_layer1, " direct targets | ",
          "Secondary TFs: ", ns$n_secondary, "\n",
          "Outer ring: ", ns$n_layer2, " downstream targets"
        ),
        caption = "Node shapes: diamond = TF, circle = gene | Size indicates layer"
      )
  } else {
    p <- p +
      labs(
        title = paste(primary_tf, "Regulatory Cascade"),
        subtitle = paste0(
          "Contrast: ", contrast_label %||% "Contrast", "\n",
          "Layer 1: ", ns$n_layer1, " direct targets | ",
          "Secondary TFs: ", ns$n_secondary,
          ifelse(ns$n_secondary > 0,
                 paste0(" (", paste(ns$secondary_tfs, collapse = ", "), ")"), ""),
          "\n",
          "Layer 2: ", ns$n_layer2, " downstream targets"
        ),
        caption = paste0(
          "Node shapes: ", intToUtf8(9670), " = TF, ", intToUtf8(9679), " = Gene\n",
          "Node size indicates regulatory layer (larger = closer to ", primary_tf, ")\n",
          "Edge style: dotted = weak binding, solid = strong binding\n",
          "Edge width: thicker = greater chromatin accessibility change"
        )
      )
  }

  # --- Save ---
  suffix <- paste0(
    ifelse(circular, "_circular", ""),
    ifelse(compact, "_compact", "")
  )
  save_cascade_plot(p, output_prefix, suffix = suffix,
                    save_pdf = save_pdf, save_png = save_png,
                    width = width, height = height, dpi = dpi)

  return(p)
}

# ==============================================================================
# PUBLIC API: Thin wrappers preserving original function signatures
# ==============================================================================

plot_tf_cascade <- function(g, primary_tf, contrast_label = NULL,
                            layout_type = "sugiyama",
                            output_prefix = NULL,
                            save_pdf = TRUE,
                            save_png = TRUE,
                            width = 14,
                            height = 12,
                            dpi = 300) {
  #' Create cascade network visualization (non-circular, full size)
  #'
  #' @param g igraph object from build_cascade_graph
  #' @param primary_tf Name of primary TF (center of cascade)
  #' @param contrast_label Label for contrast (shown in subtitle)
  #' @param layout_type Layout algorithm: "sugiyama", "tree", "fr", "kk"
  #' @param output_prefix File prefix for saving (without extension)

  plot_tf_cascade_core(
    g = g, primary_tf = primary_tf, contrast_label = contrast_label,
    layout_type = layout_type, output_prefix = output_prefix,
    save_pdf = save_pdf, save_png = save_png,
    width = width, height = height, dpi = dpi,
    circular = FALSE, compact = FALSE
  )
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

  plot_tf_cascade_core(
    g = g, primary_tf = primary_tf, contrast_label = contrast_label,
    output_prefix = output_prefix,
    save_pdf = save_pdf, save_png = save_png,
    width = width, height = height, dpi = dpi,
    circular = TRUE, compact = FALSE
  )
}

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

  plot_tf_cascade_core(
    g = g, primary_tf = primary_tf, contrast_label = contrast_label,
    output_prefix = output_prefix,
    save_pdf = save_pdf, save_png = save_png,
    width = width, height = height, dpi = dpi,
    circular = FALSE, compact = TRUE
  )
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

  plot_tf_cascade_core(
    g = g, primary_tf = primary_tf, contrast_label = contrast_label,
    output_prefix = output_prefix,
    save_pdf = save_pdf, save_png = save_png,
    width = width, height = height, dpi = dpi,
    circular = TRUE, compact = TRUE
  )
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
