library(tidyverse)
library(rhdf5)
library(here)

pca_path <- here::here("data/interim/FeatureAnalysis/line_differences/human/all_drugs/results/ReducedFeaturesPCA_all_drugs_human.h5")
h5ls(pca_path)
hdf5_pca <- h5read(pca_path, "features_organoids") %>% as_tibble()
hdf5_pca_meta <- h5read(pca_path, "metadata_organoids") %>% as.matrix %>% t() %>% as_tibble() %>% magrittr::set_colnames(h5read(pca_path, "metadata_names_organoids"))
hdf5_pca <- cbind(hdf5_pca_meta, hdf5_pca) %>% 
  arrange(Line, Plate, Well) %>%
  mutate(morphological_class = case_when(Line %in% c("D046T01",
                                                     "D054T01",
                                                     "D055T01") ~ "disorganized",
                                         Line %in% c("D018T01") ~ "organized",
                                         Line %in% c("D027T01",
                                                     "D013T01",
                                                     "D030T01") ~ "weakly organized",
                                         Line %in% c("D010T01",
                                                     "D007T01") ~ "weakly disorganized",
                                         Line %in% c("D004T01",
                                                     "D019T01",
                                                     "D020T01",
                                                     "D022T01") ~ "intermediate",
                                         Line %in% c("D018T01") ~ "organized",
                                         TRUE ~ "other") %>% factor(levels = c("organized", "weakly organized",
                                                                               "intermediate", "weakly disorganized",
                                                                               "disorganized", "other"))) %>%
  janitor::clean_names() %>%
  magrittr::set_colnames(colnames(.) %>% str_replace('v', 'PC'))
  
print(head(hdf5_pca))

saveRDS(hdf5_pca, here("data/interim/PhenotypeSpectrum/hdf5_pca_absolute_all_drugs.Rds"))

# 