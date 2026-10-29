#!/usr/bin/env bash
# Rebuild docs/CNV_CALLING_EXPLAINED.pdf from its .qmd source.
#
# Two steps and no new dependencies: quarto (already required by every
# notebook stage) renders the HTML, and headless Chrome prints it. A LaTeX
# toolchain would be a ~100 MB download for a document that is mostly tables,
# and quarto's own PDF path needs one.
#
# The intermediate HTML is deliberately NOT committed -- it is a 170 KB
# self-contained build artefact of a 550-line source file.
set -euo pipefail

cd "$(dirname "$0")"
SRC=CNV_CALLING_EXPLAINED
CHROME="${CHROME:-/Applications/Google Chrome.app/Contents/MacOS/Google Chrome}"

command -v quarto >/dev/null || { echo "quarto not on PATH" >&2; exit 1; }
[ -x "$CHROME" ] || { echo "Chrome not found at: $CHROME (set \$CHROME)" >&2; exit 1; }

quarto render "${SRC}.qmd" --to html

# --virtual-time-budget lets quarto's own JS settle before the snapshot; without
# it the table of contents is sometimes printed mid-relocation.
"$CHROME" --headless --disable-gpu --no-pdf-header-footer \
          --virtual-time-budget=10000 \
          --print-to-pdf="${SRC}.pdf" "file://$PWD/${SRC}.html" 2>/dev/null

rm -f "${SRC}.html"
echo "wrote $PWD/${SRC}.pdf"
