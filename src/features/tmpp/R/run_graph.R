library(tidyverse)
library(here)
#library(feather)
library(monocle3)
#library(ggrastr)
library(microbenchmark)

# loading data
ods <- readRDS(here("R/ods_clust_50"))
print("loaded data")

# defining functions
iterate_length <- function(length_in){
  print(paste0("k is ", length_in))
  ods_graph = learn_graph(ods, verbose = TRUE, 
                              close_loop = TRUE, 
                              learn_graph_control = list(minimal_branch_len = length_in))
  
  saveRDS(ods_graph, paste0("ods_graph_", length_in))
}

diagnose_l_iter <- function(length_in){
  print(paste0("diagnosing length ", length_in))
  bm <- microbenchmark(iterate_length(length_in), times = 1)
  return(bm)
}

# running benchmark
bm_l <- lapply(c(10, 20, 30, 40, 50, 100, 1000, 10000), diagnose_l_iter)
saveRDS(bm_l, "benchmark_graph_run.Rds")