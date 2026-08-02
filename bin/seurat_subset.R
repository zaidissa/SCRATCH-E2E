#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(readr)
library(tibble)
library(optparse)

#Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")

#

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "Input dump name", metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character",
              help = "Metadata with cell status (TME/Malignant)", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = './',
              help = "Output directory path", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#

if (length(opt) == 3) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call. = FALSE)
}

#

seurat_object <- readRDS(file = opt$file)
metadata <- readr::read_csv(opt$metadata)

metadata <- metadata %>%
  tibble::column_to_rownames(var = "barcode") %>%
  as.data.frame()

seurat_object <- Seurat::AddMetaData(
  seurat_object, metadata = metadata
  )

seurat_object <- subset(
  seurat_object, subset = cell_status == "TME"
)

#

saveRDS(
  seurat_object, file = file.path(opt$outdir, gsub(".RDS", "_filtered.RDS", basename(opt$file))))
