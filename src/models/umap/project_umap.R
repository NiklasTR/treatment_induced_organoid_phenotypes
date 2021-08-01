# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==3) {
  stop("Arguments must be supplied (raw input file name, harmony input file name, seed).n", call.=FALSE)
} 

print(args)

# seed
set.seed(args[3])

# packages
library(harmony)
library(tidyr)
library(tibble)
library(dplyr)
library(readr)
library(stringr)
library(magrittr)
library(monocle3)
library(readxl)
print("loaded libraries")

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==3) {
  stop("Arguments must be supplied (filename harmony, filename raw, seed).n", call.=FALSE)
} 

PATH = paste0(here::here(), "/")

# read object from file
obj <- readRDS(paste0(PATH, args[1]))
obj_h <- readRDS(paste0(PATH, args[2]))

## Run UMAP embedding
print("starting UMAP embedding")
# harmony
obj <- reduce_dimension(obj,
                        reduction_method = "UMAP",
                        umap.min_dist = 0.1,
                        umap.n_neighbors = 15L,
                        umap.fast_sgd=TRUE,
                        cores = parallel::detectCores(),
                        verbose = TRUE)
# non-harmony
obj_h <- reduce_dimension(obj_h,
                        reduction_method = "UMAP",
                        umap.min_dist = 0.1,
                        umap.n_neighbors = 15L,
                        umap.fast_sgd=TRUE, 
                        cores = parallel::detectCores(),
                        verbose = TRUE)

# Save harmony result
saveRDS(obj, args[1])
saveRDS(obj_h, args[2])
print("saved UMAP projected objects at:")
print(args[1])
print(args[2])
