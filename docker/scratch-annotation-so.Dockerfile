# ---------------------------------------------------------------------------
# scratch-annotation, with SeuratObject pinned — the fix for the pairing that
# has been breaking the annotation stage since the 2026-02-13 `latest` build.
#
# `latest` ships Seurat 5.0.0 against SeuratObject 5.3.0. Seurat 5.0.0 still
# calls `PackageCheck()` and `FetchData(slot=)` internally, and SeuratObject
# made both DEFUNCT after 5.2.0, so the stage dies inside Seurat's own code:
#
#     ANNOTATION:SCEASY_CONVERTER_ONE / SCYTPE_MAJOR_ANNOTATION
#     Error: ! `PackageCheck()` was deprecated in SeuratObject 5.0.0 and is
#            now defunct.
#       1. Seurat::as.SingleCellExperiment(seurat_object)
#       3.   SeuratObject::PackageCheck("SingleCellExperiment", error = FALSE)
#
# (Cirro, 2026-10-16: five failed tasks, three retries each, run aborted.)
# No notebook can work around it — nothing in this repository calls either
# function; the skew is between two packages inside the image.
#
# 5.2.0 and not the newest, and not an upgrade of Seurat: it is the version
# syedsazaidi/scratch-cluster:latest carries against the same Seurat 5.0.0, and
# that pairing runs the clustering notebook (FeaturePlot, VlnPlot, FindMarkers)
# in production. Pinning to it makes the two images agree instead of inventing
# a third combination. docker/scratch-annotation.Dockerfile asks for Seurat
# v5.3.1 but the published image reports 5.0.0, so a from-scratch rebuild would
# have to be re-verified against every annotation notebook; this does not.
#
#   docker build --platform=linux/amd64 \
#       -f docker/scratch-annotation-so.Dockerfile \
#       -t syedsazaidi/scratch-annotation:seuratobject-5.2.0 .
#   docker push syedsazaidi/scratch-annotation:seuratobject-5.2.0
# ---------------------------------------------------------------------------
FROM --platform=linux/x86_64 syedsazaidi/scratch-annotation:latest

# dependencies = FALSE: everything SeuratObject 5.2.0 needs (Matrix, sp, rlang,
# spam) is already in the image at a satisfying version, and letting the install
# upgrade them is how the skew got here in the first place.
RUN Rscript -e "install.packages( \
      'https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_5.2.0.tar.gz', \
      repos = NULL, type = 'source', dependencies = FALSE)"

# Fail the BUILD, not a pipeline task half an hour in. These are the exact
# calls that aborted on Cirro, plus the two plotting functions named in the
# warning at the top of conf/containers.config.
RUN Rscript -e "suppressPackageStartupMessages({ \
      library(Seurat); library(SeuratObject); library(SingleCellExperiment); \
      library(sceasy); library(reticulate) }); \
      stopifnot(packageVersion('SeuratObject') == '5.2.0'); \
      invisible(SeuratObject::PackageCheck('SingleCellExperiment', error = FALSE)); \
      set.seed(1); \
      m <- matrix(rpois(4000, 3), nrow = 100, \
                  dimnames = list(paste0('G', 1:100), paste0('C', 1:40))); \
      so <- CreateSeuratObject(counts = m); \
      so <- NormalizeData(so, verbose = FALSE); \
      so <- FindVariableFeatures(so, nfeatures = 50, verbose = FALSE); \
      so <- ScaleData(so, verbose = FALSE); \
      so <- RunPCA(so, npcs = 5, verbose = FALSE); \
      so <- FindNeighbors(so, dims = 1:5, verbose = FALSE); \
      so <- FindClusters(so, verbose = FALSE); \
      invisible(FetchData(so, vars = 'G1', slot = 'data')); \
      print(class(FeaturePlot(so, features = 'G1'))); \
      print(class(VlnPlot(so, features = 'G1'))); \
      sce <- as.SingleCellExperiment(so); \
      use_python('/opt/venv/bin/python', required = TRUE); \
      convertFormat(sce, from = 'sce', to = 'anndata', outFile = '/tmp/check.h5ad'); \
      stopifnot(file.exists('/tmp/check.h5ad'), file.size('/tmp/check.h5ad') > 0); \
      cat('OK  Seurat', as.character(packageVersion('Seurat')), \
          '| SeuratObject', as.character(packageVersion('SeuratObject')), \
          '| h5ad', file.size('/tmp/check.h5ad'), 'bytes\n')" \
  && rm -f /tmp/check.h5ad

CMD ["bash"]
