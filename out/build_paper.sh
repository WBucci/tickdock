#!/bin/bash
# Rebuild the submission document. Requires pandoc + LibreOffice on Windows.
cd "$(dirname "$0")/.."
PANDOC="${PANDOC:-$(command -v pandoc || echo "$LOCALAPPDATA/Pandoc/pandoc.exe")}"
SOFFICE="/c/Program Files/LibreOffice/program/soffice.exe"
# figures: wsl python3 scripts... (see out/make_figures.py)
"$PANDOC" docs/paper_A_submission.md -o out/paper_A.docx --resource-path=.
"$SOFFICE" --headless --convert-to pdf --outdir out out/paper_A.docx
