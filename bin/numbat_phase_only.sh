#!/usr/bin/env bash
# numbat_phase_only.sh <pileup_and_phase.R> <out.R>
#
# Writes a copy of Numbat's own pileup_and_phase.R with its cellsnp-lite step
# removed, so the joint genotyping, Eagle2 phasing and allele-table steps can
# run on pileups that were made earlier, one library per task.
#
# Why derive it instead of keeping a copy in bin/: everything after the pileup
# (numbat:::genotype, the Eagle2 commands, numbat:::preprocess_allele, the
# 0|1 / 1|0 filter) stays exactly the code of the Numbat version in the image.
# A vendored copy would drift away from that version without telling anyone.
#
# What is cut is one block, from
#     cat('Running pileup\n')
# through
#     system(glue('sh {script} 2>&1 | tee {outdir}/pileup.log'), intern = FALSE)
# Those lines only write run_pileup.sh and run it. The lines above them only
# build the command strings (in 10x mode they touch no files), and everything
# below them reads {outdir}/pileup/{sample}/cellSNP.*. If either marker is not
# found exactly once, the upstream script has changed shape and this stops,
# rather than run a script that either piles up again or skips something else.
set -euo pipefail

src="$1"; out="$2"
[ -s "$src" ] || { echo "numbat_phase_only: '$src' is missing or empty" >&2; exit 1; }

awk '
    index($0, "cat(\x27Running pileup") == 1                 { skip = 1; start++ }
    !skip                                                    { print }
    skip && index($0, "system(glue(\x27sh {script} 2>&1 | tee {outdir}/pileup.log") == 1 { skip = 0; stop++ }
    END {
        if (start != 1 || stop != 1) {
            printf "numbat_phase_only: expected one pileup block, found start=%d stop=%d\n", start, stop > "/dev/stderr"
            exit 3
        }
    }
' "$src" > "$out"

# Check the result as well as the markers: no pileup may run, and the steps
# that have to run must still be there.
if grep -q "pileup.log" "$out"; then
    echo "numbat_phase_only: the cellsnp-lite call is still in $out" >&2; exit 3
fi
for keep in "numbat:::genotype" "run_phasing.sh" "numbat:::preprocess_allele"; do
    grep -qF "$keep" "$out" || {
        echo "numbat_phase_only: '$keep' is missing from $out" >&2; exit 3; }
done
