#!/bin/bash
################################################################################
# RNA-ATAC-TF-Network Pipeline Runner
################################################################################
#
# A multi-omics integration tool that combines TOBIAS TF footprinting,
# ATAC-seq accessibility, and RNA-seq expression to build pathway-centric
# transcription factor regulatory networks.
#
# Usage:
#   ./run_pipeline.sh <contrast> [options]
#
# Examples:
#   ./run_pipeline.sh NvT              # Run for Normal vs Transformed contrast
#   ./run_pipeline.sh CvNC             # Run for Canonical vs NonCanonical
#   ./run_pipeline.sh all              # Run all contrasts defined in config
#   ./run_pipeline.sh tf_focus         # Run TF-focused network only
#   ./run_pipeline.sh NvT --tf_focus   # Run NvT + TF-focused analysis
#   ./run_pipeline.sh --help           # Show this help message
#
# For more information, see: docs/WORKFLOW.md
#
################################################################################

VERSION="1.0.0"

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG="$SCRIPT_DIR/config/config.yaml"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo ""
    echo -e "${BLUE}============================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================================${NC}"
}

print_step() {
    echo -e "${GREEN}>>> $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}WARNING: $1${NC}"
}

print_error() {
    echo -e "${RED}ERROR: $1${NC}"
}

show_help() {
    cat << EOF
RNA-ATAC-TF-Network Pipeline v${VERSION}

DESCRIPTION:
  Multi-omics integration tool that combines TOBIAS TF footprinting,
  ATAC-seq accessibility, and RNA-seq expression to build pathway-centric
  transcription factor regulatory networks.

USAGE:
  ./run_pipeline.sh <contrast|command> [options]

COMMANDS:
  <contrast>      Run pipeline for specified contrast (e.g., NvT, CvNC)
  all             Run pipeline for all contrasts defined in config
  tf_focus        Run TF-focused network analysis only

OPTIONS:
  --tf_focus      Also run TF-focused analysis after main pipeline
  --config FILE   Use alternate config file (default: config/config.yaml)
  --help, -h      Show this help message
  --version, -v   Show version information

EXAMPLES:
  ./run_pipeline.sh NvT                  # Run Normal vs Transformed
  ./run_pipeline.sh CvNC                 # Run Canonical vs NonCanonical
  ./run_pipeline.sh all                  # Run all contrasts
  ./run_pipeline.sh tf_focus             # TF-focused networks only
  ./run_pipeline.sh NvT --tf_focus       # NvT + TF-focused
  ./run_pipeline.sh NvT --config my.yaml # Use custom config

OUTPUT:
  Tables:    ./output/<contrast>/tables/
  Plots:     ./output/<contrast>/plots/
  Cytoscape: ./output/<contrast>/cytoscape/
  TF-Focus:  ./output/tf_focus/

For detailed documentation, see:
  - docs/WORKFLOW.md         Step-by-step usage guide
  - docs/INPUT_DATA_FORMAT.md  Input file specifications
  - docs/OUTPUT_FILES.md     Output file descriptions

EOF
}

show_version() {
    echo "RNA-ATAC-TF-Network Pipeline v${VERSION}"
}

# Check dependencies
check_dependencies() {
    print_header "Checking dependencies..."

    # Check Python
    if ! command -v python3 &> /dev/null; then
        print_error "Python 3 not found. Please install Python 3."
        echo "  Run: ./install/setup.sh"
        exit 1
    fi
    echo "  Python: $(python3 --version)"

    # Check required Python packages
    python3 -c "import pandas, numpy, yaml, scipy" 2>/dev/null || {
        print_warning "Some Python packages may be missing."
        echo "  Run: pip install -r install/requirements.txt"
        exit 1
    }

    # Check R
    if ! command -v Rscript &> /dev/null; then
        print_error "R not found. Please install R."
        echo "  Run: ./install/setup.sh"
        exit 1
    fi
    echo "  R: $(Rscript --version 2>&1 | head -1)"

    # Check required R packages
    Rscript -e "suppressPackageStartupMessages({
        required <- c('argparse', 'yaml', 'dplyr', 'tidyr', 'ggplot2',
                      'ggraph', 'igraph', 'msigdbr', 'stringr', 'viridis')
        missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
        if (length(missing) > 0) {
            cat('Missing R packages:', paste(missing, collapse=', '), '\n')
            quit(status = 1)
        }
    })" 2>/dev/null || {
        print_warning "Some R packages may be missing."
        echo "  Run: Rscript install/install_r_packages.R"
        exit 1
    }

    echo -e "  ${GREEN}All dependencies OK!${NC}"
}

# Check config file exists
check_config() {
    if [ ! -f "$CONFIG" ]; then
        print_error "Config file not found: $CONFIG"
        echo ""
        echo "  Please create a config file:"
        echo "    1. Copy config/config_template.yaml to config/config.yaml"
        echo "    2. Update paths and settings for your data"
        echo ""
        exit 1
    fi
}

# Run pipeline for a single contrast
run_contrast() {
    local CONTRAST=$1

    print_header "Running pipeline for contrast: $CONTRAST"

    # Step 1: Build TF-gene network
    print_step "Step 1: Building TF-gene edges..."
    python3 "$SCRIPT_DIR/scripts/build_tf_gene_network.py" \
        --config "$CONFIG" \
        --contrast "$CONTRAST"

    # Find output directory
    local OUT_DIR="$SCRIPT_DIR/output/$CONTRAST"
    local EDGES_FILE="$OUT_DIR/tables/tf_gene_edges_${CONTRAST}.tsv"
    local NODES_FILE="$OUT_DIR/tables/tf_gene_nodes_${CONTRAST}.tsv"

    if [ ! -f "$EDGES_FILE" ]; then
        print_error "Edge file not found: $EDGES_FILE"
        exit 1
    fi

    # Step 2: Generate network plots
    print_step "Step 2: Generating pathway network plots..."
    Rscript "$SCRIPT_DIR/scripts/plot_pathway_network.R" \
        --edges "$EDGES_FILE" \
        --nodes "$NODES_FILE" \
        --contrast "$CONTRAST" \
        --config "$CONFIG" \
        --all_pathways \
        --emap

    print_step "Contrast $CONTRAST complete!"
    echo "  Outputs saved to: $OUT_DIR"
}

# Run TF-focused network analysis
run_tf_focus() {
    local CONTRAST=${1:-NvT}

    print_header "Running TF-Focused Network Analysis"

    # Check if TF focus is enabled in config
    TF_FOCUS_ENABLED=$(python3 -c "import yaml; cfg=yaml.safe_load(open('$CONFIG')); print(cfg.get('tf_focus', {}).get('enabled', False))")

    if [ "$TF_FOCUS_ENABLED" != "True" ]; then
        print_warning "TF focus mode is disabled in config.yaml"
        print_step "To enable, set tf_focus.enabled: true in config.yaml"
        return 0
    fi

    # Check if prerequisite files exist
    local RNA_DE_FILE="$SCRIPT_DIR/output/$CONTRAST/tables/rna_de_${CONTRAST}.tsv"
    local PEAK2GENE_FILE="$SCRIPT_DIR/output/$CONTRAST/tables/peak2gene_${CONTRAST}.tsv"

    if [ ! -f "$RNA_DE_FILE" ] || [ ! -f "$PEAK2GENE_FILE" ]; then
        print_warning "Required files not found. Running base pipeline first..."
        run_contrast "$CONTRAST"
    fi

    # Step 1: Build TF-focused network
    print_step "Step 1: Building TF-focused cascade network..."
    python3 "$SCRIPT_DIR/scripts/build_tf_focus_network.py" \
        --config "$CONFIG" \
        --contrast "$CONTRAST"

    # Find output files
    local TF_FOCUS_DIR="$SCRIPT_DIR/output/tf_focus"

    if [ ! -d "$TF_FOCUS_DIR/tables" ]; then
        print_error "TF focus tables not generated"
        return 1
    fi

    # Step 2: Generate visualizations for each TF
    print_step "Step 2: Generating TF-focused network plots..."

    # Get TF names from config
    TF_NAMES=$(python3 -c "
import yaml
cfg = yaml.safe_load(open('$CONFIG'))
tfs = cfg.get('tf_focus', {}).get('tfs', [])
for tf in tfs:
    print(tf.get('name', ''))
")

    # Get contrast label from config
    CONTRAST_LABEL=$(python3 -c "
import yaml
cfg = yaml.safe_load(open('$CONFIG'))
print(cfg.get('tf_focus', {}).get('contrast_label', 'Contrast'))
")

    for TF_NAME in $TF_NAMES; do
        if [ -z "$TF_NAME" ]; then
            continue
        fi

        # Find the edges file for this TF
        EDGES_PATTERN="$TF_FOCUS_DIR/tables/tf_focus_edges_${TF_NAME}_*.tsv"
        EDGES_FILE=$(ls $EDGES_PATTERN 2>/dev/null | head -1)

        if [ -z "$EDGES_FILE" ] || [ ! -f "$EDGES_FILE" ]; then
            print_warning "No edges file found for $TF_NAME"
            continue
        fi

        # Get matching nodes file
        NODES_FILE="${EDGES_FILE/edges/nodes}"

        if [ ! -f "$NODES_FILE" ]; then
            print_warning "No nodes file found for $TF_NAME"
            continue
        fi

        print_step "  Plotting $TF_NAME cascade network..."

        Rscript "$SCRIPT_DIR/scripts/plot_tf_focus_network.R" \
            --edges "$EDGES_FILE" \
            --nodes "$NODES_FILE" \
            --tf "$TF_NAME" \
            --contrast "$CONTRAST_LABEL" \
            --config "$CONFIG" \
            --circular \
            --compact
    done

    print_step "TF-focused analysis complete!"
    echo "  Outputs saved to: $TF_FOCUS_DIR"
}

# Main
main() {
    cd "$SCRIPT_DIR"

    # Parse arguments
    local RUN_TF_FOCUS=false
    local CONTRAST=""
    local SHOW_HELP=false
    local SHOW_VERSION=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            --help|-h)
                SHOW_HELP=true
                shift
                ;;
            --version|-v)
                SHOW_VERSION=true
                shift
                ;;
            --tf_focus)
                RUN_TF_FOCUS=true
                shift
                ;;
            --config)
                CONFIG="$2"
                shift 2
                ;;
            tf_focus)
                RUN_TF_FOCUS=true
                shift
                ;;
            all|NvT|CvNC|*)
                if [[ "$1" != --* ]]; then
                    CONTRAST="$1"
                fi
                shift
                ;;
        esac
    done

    # Handle help and version
    if [ "$SHOW_HELP" = true ]; then
        show_help
        exit 0
    fi

    if [ "$SHOW_VERSION" = true ]; then
        show_version
        exit 0
    fi

    # Show usage if no arguments
    if [ -z "$CONTRAST" ] && [ "$RUN_TF_FOCUS" = false ]; then
        show_help
        exit 0
    fi

    # Run pipeline
    check_config
    check_dependencies

    # Run main pipeline
    if [ "$CONTRAST" == "all" ]; then
        # Get all contrasts from config
        CONTRASTS=$(python3 -c "
import yaml
cfg = yaml.safe_load(open('$CONFIG'))
for c in cfg.get('contrasts', {}).keys():
    print(c)
")
        for c in $CONTRASTS; do
            run_contrast "$c"
        done
    elif [ -n "$CONTRAST" ] && [ "$CONTRAST" != "tf_focus" ]; then
        run_contrast "$CONTRAST"
    fi

    # Run TF-focused analysis
    if [ "$RUN_TF_FOCUS" == "true" ]; then
        local TF_CONTRAST=${CONTRAST:-NvT}
        if [ "$TF_CONTRAST" == "all" ]; then
            TF_CONTRAST="NvT"
        fi
        run_tf_focus "$TF_CONTRAST"
    fi

    print_header "Pipeline complete!"
    echo ""
    echo "Output files:"
    if [ -n "$CONTRAST" ] && [ "$CONTRAST" != "tf_focus" ]; then
        echo "  Tables:    ./output/$CONTRAST/tables/"
        echo "  Plots:     ./output/$CONTRAST/plots/"
        echo "  Cytoscape: ./output/$CONTRAST/cytoscape/"
    fi
    if [ "$RUN_TF_FOCUS" == "true" ]; then
        echo "  TF-Focus:  ./output/tf_focus/"
    fi
    echo ""
}

main "$@"
