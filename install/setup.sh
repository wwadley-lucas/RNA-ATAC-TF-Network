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

# Verify minimum Python version (3.8+)
python3 -c "import sys; assert sys.version_info >= (3, 8), f'Python 3.8+ required, found {sys.version}'" 2>/dev/null || {
    echo -e "${RED}ERROR: Python 3.8+ is required. Found: $PYTHON_VERSION${NC}"
    exit 1
}

# Check R
echo ""
echo "Checking R installation..."
if ! command -v Rscript &> /dev/null; then
    echo -e "${RED}ERROR: R not found. Please install R 4.0+${NC}"
    exit 1
fi
R_VERSION=$(Rscript --version 2>&1 | head -1)
echo -e "${GREEN}  Found: $R_VERSION${NC}"

# Verify minimum R version (4.0+)
Rscript -e "if (getRversion() < '4.0') { cat('ERROR: R 4.0+ required\n'); quit(status=1) }" 2>/dev/null || {
    echo -e "${RED}ERROR: R 4.0+ is required. Found: $R_VERSION${NC}"
    exit 1
}

# Check for external tools (TOBIAS, HOMER)
echo ""
echo "============================================================"
echo "Checking external tool versions..."
echo "============================================================"

# Check TOBIAS
if command -v TOBIAS &> /dev/null; then
    TOBIAS_VERSION=$(TOBIAS --version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+')
    TOBIAS_MAJOR=$(echo "$TOBIAS_VERSION" | cut -d. -f1)
    TOBIAS_MINOR=$(echo "$TOBIAS_VERSION" | cut -d. -f2)
    if [ "$TOBIAS_MAJOR" -eq 0 ] && [ "$TOBIAS_MINOR" -lt 14 ] 2>/dev/null; then
        echo -e "${YELLOW}WARNING: TOBIAS version $TOBIAS_VERSION detected. Version 0.14+ is required.${NC}"
        echo -e "${YELLOW}  BINDetect output format changed in 0.14. Please upgrade: pip install tobias>=0.14${NC}"
    else
        echo -e "${GREEN}  Found TOBIAS $TOBIAS_VERSION (>= 0.14 required)${NC}"
    fi
else
    echo -e "${YELLOW}WARNING: TOBIAS not found. Install with: pip install tobias>=0.14${NC}"
    echo -e "${YELLOW}  TOBIAS is required for TF footprinting (ATACorrect, FootprintScores, BINDetect)${NC}"
fi

# Check HOMER
if command -v annotatePeaks.pl &> /dev/null; then
    HOMER_VERSION=$(annotatePeaks.pl 2>&1 | head -1 | grep -oE 'v[0-9]+\.[0-9]+' || echo "unknown")
    echo -e "${GREEN}  Found HOMER $HOMER_VERSION (v4.11+ recommended)${NC}"
else
    echo -e "${YELLOW}WARNING: HOMER (annotatePeaks.pl) not found.${NC}"
    echo -e "${YELLOW}  HOMER is required for peak annotation. Install from: http://homer.ucsd.edu/homer/introduction/install.html${NC}"
    echo -e "${YELLOW}  Alternatively, ChIPseeker can be used if output format matches.${NC}"
fi

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
