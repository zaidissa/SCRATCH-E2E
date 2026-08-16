#!/usr/bin/env bash
# Fetch the three external reference files Numbat needs.
#
#   bin/fetch_numbat_refs.sh /path/to/refs [hg38|hg19]
#
# Downloads ~10 GB. Run once to a shared location; on HPC put it somewhere
# group-readable so the cohort does not each keep a copy.
set -euo pipefail

DEST="${1:?usage: fetch_numbat_refs.sh <dest_dir> [hg38|hg19]}"
BUILD="${2:-hg38}"
case "$BUILD" in hg38|hg19) ;; *) echo "build must be hg38 or hg19" >&2; exit 2 ;; esac
mkdir -p "$DEST"; cd "$DEST"

echo "==> 1/3 SNP VCF ($BUILD)"
VCF="genome1K.phase3.SNP_AF5e2.chr1toX.${BUILD}.vcf.gz"
[ -f "$VCF" ] || curl -fL --retry 3 -o "$VCF" \
  "https://sourceforge.net/projects/cellsnp/files/SNPlist/${VCF}/download"

echo "==> 2/3 1000G phasing panel ($BUILD)"
if [ ! -d "1000G_${BUILD}" ]; then
  curl -fL --retry 3 -o "1000G_${BUILD}.zip" \
    "http://pklab.med.harvard.edu/teng/data/1000G_${BUILD}.zip"
  unzip -q "1000G_${BUILD}.zip" && rm -f "1000G_${BUILD}.zip"
fi

echo "==> 3/3 Eagle2 genetic map"
# The map ships inside the Eagle2 tarball; pull just the tables directory.
if [ ! -f "genetic_map_${BUILD}_withX.txt.gz" ]; then
  curl -fL --retry 3 -o Eagle.tar.gz \
    "https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz"
  tar -xzf Eagle.tar.gz --strip-components=2 \
      "Eagle_v2.4.1/tables/genetic_map_${BUILD}_withX.txt.gz"
  rm -f Eagle.tar.gz
fi

cat <<MSG

Done. Run with:

  nextflow run . -profile docker \\
    --run_numbat true \\
    --numbat_gmap     $DEST/genetic_map_${BUILD}_withX.txt.gz \\
    --numbat_snpvcf   $DEST/$VCF \\
    --numbat_paneldir $DEST/1000G_${BUILD} \\
    --numbat_genome   ${BUILD}

The build must match --genome. An hg19 panel against hg38 BAMs does not error --
it silently produces wrong calls.
MSG
