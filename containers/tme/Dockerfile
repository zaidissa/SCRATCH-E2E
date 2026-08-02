# SCRATCH-E2E — TME characterisation image (default).
#
# Adds the R packages the TME notebook needs for compositional differential
# abundance to the shared analysis image. propeller (speckle) is the DEFAULT
# compositional test: it stays in R, runs inside the notebook, and adds a few MB
# rather than the ~2 GB that scCODA's TensorFlow dependency costs.
#
# scCODA remains available as an opt-in confirmatory test -- see
# containers/tme/Dockerfile.scoda, which builds ON TOP of this image.
#
# Build (base is linux/amd64; state it explicitly so this reproduces on an
# Apple-silicon workstation):
#
#   docker build --platform linux/amd64 -f containers/tme/Dockerfile \
#     -t syedsazaidi/scratch-tme:latest .

FROM syedsazaidi/scratch-cluster:latest

LABEL org.opencontainers.image.title="scratch-tme" \
      org.opencontainers.image.description="SCRATCH-E2E TME characterisation: analysis image plus speckle/propeller for compositional differential abundance." \
      org.opencontainers.image.source="https://github.com/syedsazaidi/SCRATCH-E2E"

# speckle needs edgeR, which the base image does not carry. limma and statmod
# are already present.
RUN set -eux; \
    Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) \
                  install.packages("BiocManager", repos="https://cloud.r-project.org"); \
                BiocManager::install(c("edgeR","speckle"), ask=FALSE, update=FALSE)'

# Unlike the TensorFlow check in Dockerfile.scoda, this one runs natively under
# emulation, so it is a hard gate: a missing package here means the notebook
# silently falls back to reporting no compositional test at all.
RUN set -eux; \
    Rscript -e 'stopifnot(requireNamespace("speckle", quietly=TRUE), \
                          requireNamespace("edgeR",  quietly=TRUE)); \
                cat("VERIFIED speckle", as.character(packageVersion("speckle")), \
                    "| edgeR", as.character(packageVersion("edgeR")), "\n")'
