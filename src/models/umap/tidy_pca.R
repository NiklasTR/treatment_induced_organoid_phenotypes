# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

library(tidyverse)
library(rhdf5)
library(here)

pca_path <- here::here("data/interim/FeatureAnalysis/line_differences/human/all_drugs/results/ReducedFeaturesPCA_all_drugs_human.h5")
h5ls(pca_path)
hdf5_pca <- h5read(pca_path, "features_organoids") %>% as_tibble()
hdf5_pca_meta <- h5read(pca_path, "metadata_organoids") %>% as.matrix %>% t() %>% as_tibble() %>% magrittr::set_colnames(h5read(pca_path, "metadata_names_organoids"))

print(hdf5_pca_meta$Line %>% table())

hdf5_pca <- cbind(hdf5_pca_meta, hdf5_pca) %>% 
  arrange(Line, Plate, Well) %>%
  janitor::clean_names() %>%
  magrittr::set_colnames(colnames(.) %>% str_replace('v', 'PC'))
  
print(head(hdf5_pca))

saveRDS(hdf5_pca, here("data/interim/PhenotypeSpectrum/hdf5_pca_absolute_all_drugs.Rds"))
