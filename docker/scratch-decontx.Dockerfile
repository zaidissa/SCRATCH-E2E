# ---------------------------------------------------------------------------
# scratch-decontx — DecontX ambient-RNA removal (--ambient_method decontx)
#
# Built ON TOP of the QC image rather than from a Bioconductor base: that image
# already carries Seurat, DropletUtils (which reads and writes cellranger h5)
# and hdf5r, which is everything bin/decontx_run.R needs apart from celda
# itself. Layering keeps the two images on identical Seurat/Matrix versions —
# a separate base would reintroduce exactly the version skew that broke
# scratch-annotation (Seurat 5.0.0 against SeuratObject 5.3.0).
#
#   docker build --platform=linux/amd64 \
#       -f docker/scratch-decontx.Dockerfile -t syedsazaidi/scratch-decontx:latest .
#   docker push syedsazaidi/scratch-decontx:latest
# ---------------------------------------------------------------------------
FROM --platform=linux/x86_64 syedsazaidi/scratch-qc:latest

RUN Rscript -e "BiocManager::install('celda', ask = FALSE, update = FALSE)"

# Fail the BUILD, not a four-hour pipeline task, if anything the script needs is
# missing from the image.
RUN Rscript -e "suppressPackageStartupMessages({ \
      library(celda); library(DropletUtils); library(SingleCellExperiment); \
      library(Matrix); library(Seurat) }); \
      stopifnot(is.function(celda::decontX), is.function(DropletUtils::write10xCounts)); \
      cat('celda', as.character(packageVersion('celda')), \
          '| Seurat', as.character(packageVersion('Seurat')), '\n')"

CMD ["bash"]
