#!/usr/bin/Rscript

# packages
library(tidyverse)
library(monocle3)

# run:
args = commandArgs(trailingOnly=TRUE)

# not run: 
# args <- c("data/interim/FeatureAnalysis/organoid_viability/human/results", "data/interim/PhenotypeSpectrum/harmony_umap_absolute_all_drugs.Rds")

print(args)

if (!length(args)==2) {
  stop("Arguments must be supplied (input directory logs, input monocle object).n", call.=FALSE)
} 



# find classification_logs and aggregate results for UMAP filtering
log <- list.files(args[1], full.names = TRUE, pattern = "log.csv") %>%
  as_tibble() %>% janitor::clean_names() %>%
  # reading logs
  mutate(data = purrr::map(value, ~ .x %>% read_csv)) %>% 
  # crunching filenames to extract organoid line name in a hacky way
  separate(value, into = c("path", "ending"), sep = "_", fill = "right") %>% 
  mutate(line = substr(ending, nchar(ending)-6, nchar(ending))) %>%
  dplyr::select(everything(), -path, -ending) %>% 
  # counting the number of objects
  mutate(n_log = purrr::map(data, ~ .x %>% nrow())) %>% 
  unnest(n_log)

# accessing the monocle object
obj <- readRDS(args[2])
umap_tidy <- reducedDims(obj)$UMAP %>% cbind(colData(obj)) %>% as_tibble() %>% janitor::clean_names()

# TEST: the number of objects in the log file has to match the number of objects in the monocle object.
obj_n <- umap_tidy %>% dplyr::count(line)
ldc_n <- log %>% dplyr::select(line, n = n_log)
stopifnot(obj_n == ldc_n)
  
# write object with added LDC information.
df <- log %>% dplyr::select(everything(), -n_log) %>% 
  unnest(data) %>% 
  dplyr::select(-line) %>% janitor::clean_names() %>%
  cbind(colData(obj), .)

# I am saving my final result
pData(obj) <- df
obj %>% saveRDS(args[2])
