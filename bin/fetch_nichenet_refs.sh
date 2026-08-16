#!/usr/bin/env bash
# Fetch NicheNet's reference networks (the repo ships unresolved Git LFS
# pointers, which fail at readRDS with an unhelpful message).
#
#   bin/fetch_nichenet_refs.sh /path/to/nichenet_resources
#
# Source: the NicheNet authors' Zenodo record for the v2 (2021) networks.
set -euo pipefail
DEST="${1:?usage: fetch_nichenet_refs.sh <dest_dir>}"
mkdir -p "$DEST"; cd "$DEST"
Z=https://zenodo.org/records/7074291/files

fetch () { # remote name -> local name expected by the notebook
  [ -s "$2" ] && { echo "  have $2"; return; }
  echo "  fetching $2"; curl -fL --retry 3 -o "$2" "$Z/$1"
}
fetch lr_network_human_21122021.rds        lr_network_human.rds
fetch ligand_target_matrix_nsga2r_final.rds ligand_target_matrix.rds
fetch weighted_networks_nsga2r_final.rds    weighted_networks.rds

echo
echo "Verifying they are real RDS, not pointers:"
for f in *.rds; do
  if head -c 24 "$f" | grep -q 'git-lfs'; then echo "  $f  STILL A POINTER"; exit 1; fi
  printf "  %-28s %s bytes\n" "$f" "$(wc -c < "$f")"
done
echo
echo "Run with:  --nichenet_assets_dir $DEST"
