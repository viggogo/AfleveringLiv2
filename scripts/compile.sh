#!/usr/bin/env bash
# Compile main.tex to PDF
# Run with: bash scripts/compile.sh

set -euo pipefail

cd "$(dirname "$0")/.."

echo "==> Removing stale auxiliary files..."
rm -f main.aux main.log main.toc main.out texput.log

echo "==> First pass..."
pdflatex -interaction=nonstopmode main.tex

echo "==> Second pass (resolves TOC & references)..."
pdflatex -interaction=nonstopmode main.tex

echo "==> Cleaning auxiliary files..."
rm -f main.aux main.log main.toc main.out texput.log

echo "==> Done! Output: main.pdf"
