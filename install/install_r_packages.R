#!/usr/bin/env Rscript
################################################################################
# RNA-ATAC-TF-Network R Package Installer
################################################################################
# Run with: Rscript install_r_packages.R
################################################################################

cat("Installing R packages for RNA-ATAC-TF-Network...\n\n")

# Required packages from CRAN
cran_packages <- c(
  "argparse",       # Command-line argument parsing
  "yaml",           # YAML configuration files
  "dplyr",          # Data manipulation
  "tidyr",          # Data tidying
  "tibble",         # Modern data frames
  "ggplot2",        # Plotting
  "ggraph",         # Network visualization
  "igraph",         # Network analysis
  "stringr",        # String manipulation
  "viridis",        # Color palettes
  "RColorBrewer",   # Color palettes
  "scales",         # Scale functions
  "tidygraph"       # Tidy network analysis
)

# Bioconductor packages
bioc_packages <- c(
  "msigdbr"         # MSigDB gene sets
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# Check for BiocManager and install if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("\nInstalling BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
}

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# Verify all packages can be loaded
cat("\n")
cat("Verifying package installation...\n")

all_packages <- c(cran_packages, bioc_packages)
failed <- character(0)

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    failed <- c(failed, pkg)
    cat(sprintf("  FAILED: %s\n", pkg))
  } else {
    cat(sprintf("  OK: %s\n", pkg))
  }
}

if (length(failed) > 0) {
  cat(sprintf("\nERROR: Failed to install %d package(s): %s\n",
              length(failed), paste(failed, collapse = ", ")))
  cat("Please install these packages manually.\n")
  quit(status = 1)
} else {
  cat("\nAll packages installed successfully!\n")
}

# Print session info
cat("\n")
cat("Session info:\n")
cat(sprintf("  R version: %s\n", R.version$version.string))
cat(sprintf("  Platform: %s\n", R.version$platform))
