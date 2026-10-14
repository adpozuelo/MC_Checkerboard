#!/usr/bin/env bash
# ==============================================================================
# Shell wrapper to reproduce all figures of the MCCB-gpu manuscript
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo "========================================================================"
echo " Reproducing All Figures for the MCCB-gpu Manuscript"
echo " Working directory: ${SCRIPT_DIR}"
echo "========================================================================"

# Check python3
if ! command -v python3 &>/dev/null; then
    echo "Error: python3 is required but not found in PATH." >&2
    exit 1
fi

python3 "${SCRIPT_DIR}/reproduce_all_figures.py" "$@"
