# pkharchenkolab/numbat-rbase + procps.
#
# Nextflow's task wrapper begins with:
#     command -v ps &>/dev/null || { >&2 echo "Command 'ps' required..."; exit 1; }
# so with trace/report/timeline enabled (they are, in nextflow.config) EVERY
# container must have `ps`. The upstream Numbat image ships no procps, so every
# NUMBAT_PILEUP task exited 1 in ~3 seconds having never run its script -- the
# only clue being a message that reads like a warning.
#
# Nothing else is changed: the Numbat R package, Eagle2, samtools, and the
# bundled hg38/hg19 references under /data all come from upstream.
FROM pkharchenkolab/numbat-rbase:latest

# --allow-releaseinfo-change: the base is Debian *testing*, and Debian has since
# promoted trixie -> forky, so a plain `apt-get update` aborts on the codename
# change (exit 100). This only accepts the rename, not new repositories.
RUN apt-get update --allow-releaseinfo-change \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*
