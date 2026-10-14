#!/usr/bin/env bash
# ==============================================================================
# Helper script to compile the MCCB-gpu Computer Physics Communications manuscript
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo "========================================================================"
echo " Compiling MCCB-gpu CPC Manuscript: mccb_gpu_cpc.tex"
echo " Directory: ${SCRIPT_DIR}"
echo "========================================================================"

# If pdflatex not in PATH, try module load if on an HPC environment with Modules
if ! command -v pdflatex &>/dev/null; then
    if command -v module &>/dev/null; then
        echo "pdflatex not found in PATH, attempting module load texlive..."
        module load texlive/20260310-GCC-13.2.0 2>/dev/null || module load texlive 2>/dev/null || true
    fi
fi

if ! command -v pdflatex &>/dev/null; then
    echo "Error: pdflatex is not installed or not in PATH." >&2
    echo "Please install TeX Live (including elsarticle and tikz) or load the appropriate module." >&2
    exit 1
fi

echo "Running pdflatex (pass 1)..."
pdflatex -interaction=nonstopmode mccb_gpu_cpc.tex > compile.log 2>&1 || true

if command -v bibtex &>/dev/null; then
    echo "Running bibtex..."
    bibtex mccb_gpu_cpc >> compile.log 2>&1 || true
fi

echo "Running pdflatex (pass 2)..."
pdflatex -interaction=nonstopmode mccb_gpu_cpc.tex >> compile.log 2>&1 || true

echo "Running pdflatex (pass 3 for references and labels)..."
pdflatex -interaction=nonstopmode mccb_gpu_cpc.tex >> compile.log 2>&1

if [ -f "mccb_gpu_cpc.pdf" ]; then
    PAGE_COUNT=$(pdfinfo mccb_gpu_cpc.pdf 2>/dev/null | grep -i Pages | awk '{print $2}' || echo "N/A")
    echo "========================================================================"
    echo " SUCCESS: mccb_gpu_cpc.pdf compiled successfully! (Pages: ${PAGE_COUNT})"
    echo "========================================================================"
else
    echo "Error: mccb_gpu_cpc.pdf was not produced. Check compile.log for errors." >&2
    exit 1
fi
