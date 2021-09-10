# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
print(.libPaths())

# libraries
library(tidyr)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(here)
#library(ggrastr)
#library(cowplot)

# read complete data
pca_tidy <- read_rds(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_tidy.Rds"))

# run queries on data
pca_tidy_aggregate <- pca_tidy %>% dplyr::filter(drug == "DMSO") %>% 
  mutate(rep = paste0("r", replicate)) %>% 
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  mutate(line = paste0(line, "_", concentration))
  #filter(size_log > 7.5) %>%
  group_by(line) %>% 
  summarise_at(vars(contains("pc")), funs(mean)) %>% 
  ungroup() %>% 
  gather(pca, value, -line) %>% rename(feature = pca, sample = line) %>% mutate(view = "morphology_view") %>% 
  mutate(feature = paste0(feature, "_", view))

# save output
pca_tidy_aggregate %>% write_rds(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_aggregate.Rds"))


