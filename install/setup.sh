#!/bin/bash
################################################################################
# RNA-ATAC-TF-Network Setup Script
################################################################################
# This script installs all required dependencies for the pipeline.
#
# Usage:
#   ./install/setup.sh
#
# Requirements:
#   - Python 3.8+
#   - R 4.0+
#   - pip
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

echo ""
echo "============================================================"
echo "RNA-ATAC-TF-Network Setup"
echo "============================================================"
echo ""

# Check Python
echo "Checking Python installation..."
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}ERROR: Python 3 not found. Please install Python 3.8+${NC}"
    exit 1
fi
PYTHON_VERSION=$(python3 --version 2>&1)
echo -e "${GREEN}  Found: $PYTHON_VERSION${NC}"

# Check R
echo ""
echo "Checking R installation..."
if ! command -v Rscript &> /dev/null; then
    echo -e "${RED}ERROR: R not found. Please install R 4.0+${NC}"
    exit 1
fi
R_VERSION=$(Rscript --version 2>&1 | head -1)
echo -e "${GREEN}  Found: $R_VERSION${NC}"

# Install Python packages
echo ""
echo "============================================================"
echo "Installing Python packages..."
echo "============================================================"
pip3 install -r "$SCRIPT_DIR/requirements.txt"

# Install R packages
echo ""
echo "============================================================"
echo "Installing R packages..."
echo "============================================================"
Rscript "$SCRIPT_DIR/install_r_packages.R"

# Make scripts executable
echo ""
echo "============================================================"
echo "Setting permissions..."
echo "============================================================"
chmod +x "$REPO_DIR/run_pipeline.sh"
chmod +x "$REPO_DIR/scripts/"*.py
chmod +x "$REPO_DIR/scripts/"*.R
echo -e "${GREEN}  Made scripts executable${NC}"

# Create output directory
mkdir -p "$REPO_DIR/output"
echo -e "${GREEN}  Created output directory${NC}"

# Summary
echo ""
echo "============================================================"
echo -e "${GREEN}Setup complete!${NC}"
echo "============================================================"
echo ""
echo "Next steps:"
echo "  1. Copy config/config_template.yaml to config/config.yaml"
echo "  2. Update config.yaml with your data paths"
echo "  3. Run the pipeline:"
echo ""
echo "     ./run_pipeline.sh YOUR_CONTRAST"
echo ""
echo "  For help:"
echo "     ./run_pipeline.sh --help"
echo ""
