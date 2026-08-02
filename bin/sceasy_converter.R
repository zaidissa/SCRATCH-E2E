#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(sceasy)
library(optparse)
library(SingleCellExperiment)
library(reticulate)
use_python("/opt/venv/bin/python", required = TRUE)
# use_python("python", required = TRUE)
# reticulate::use_python("~/scratch_pyenv/bin/python", required = TRUE)
#Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")

#

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "Input dump name", metavar = "character"),
  make_option(c("-t", "--type"), type = "character", default = "Seurat",
              help = "Dump type", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = './',
              help = "Output directory path", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
}

print(opt)

#

if(opt$type == "Seurat") {

  seurat_object <- readRDS(file = opt$file)
  sce_object <- as.SingleCellExperiment(
    seurat_object)
  
  output <- file.path(opt$outdir, gsub(".RDS", ".h5ad", basename(opt$file)))

  convertFormat(sce_object, from = "sce", to = "anndata",
                outFile = output)

}
