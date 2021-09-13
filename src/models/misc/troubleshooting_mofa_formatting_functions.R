# old code
mofa_morphology_agg <- readRDS(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_aggregate.Rds")) %>% 
  dplyr::filter(drug == "DMSO") %>% 
  mutate(rep = paste0("r", replicate)) %>% 
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  #filter(size_log > 7.5) %>%
  group_by(line) %>% 
  summarise_at(vars(contains("pc")), funs(mean)) %>% 
  ungroup() %>% 
  gather(pca, value, -line) %>% rename(feature = pca, sample = line) %>% mutate(view = "morphology_view") %>% 
  mutate(feature = paste0(feature, "_", view))

mofa_morphology_sampled <- readRDS(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_sampled.Rds")) %>% 
  dplyr::filter(drug == "DMSO") %>% 
  mutate(rep = paste0("r", replicate)) %>% 
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  #filter(size_log > 7.5) %>%
  group_by(line) %>% 
  summarise_at(vars(contains("pc")), funs(mean)) %>% 
  ungroup() %>% 
  gather(pca, value, -line) %>% rename(feature = pca, sample = line) %>% mutate(view = "morphology_view") %>% 
  mutate(feature = paste0(feature, "_", view))

# new code
mofa_morphology_agg_new <- readRDS(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_aggregate.Rds")) %>% 
  dplyr::filter(drug == "DMSO") %>% 
  mutate(rep = paste0("r", replicate)) %>% 
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  #filter(size_log > 7.5) %>%
  #group_by(line) %>% 
  #summarise_at(vars(contains("pc")), funs(mean)) %>% 
  #ungroup() %>% 
  dplyr::select(line, contains("pc")) %>%
  gather(pca, value, -line) %>% rename(feature = pca, sample = line) %>% mutate(view = "morphology_view") %>% 
  mutate(feature = paste0(feature, "_", view)) %>%
  # I average one more time over drug activity per line. This is necessary as D020 was imaged twice. In this particular case, I am creating the average activity score
  dplyr::group_by(sample, feature, view) %>% summarise(value = mean(value)) %>%
  arrange(feature, sample) %>%
  dplyr::select(sample, feature, value, view)

mofa_morphology_sampled_new <- readRDS(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_sampled.Rds")) %>% 
  dplyr::filter(drug == "DMSO") %>% 
  mutate(rep = paste0("r", replicate)) %>% 
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  #filter(size_log > 7.5) %>%
  #group_by(line) %>% 
  #summarise_at(vars(contains("pc")), funs(mean)) %>% 
  #ungroup() %>% 
  dplyr::select(line, contains("pc")) %>%
  gather(pca, value, -line) %>% rename(feature = pca, sample = line) %>% mutate(view = "morphology_view") %>% 
  mutate(feature = paste0(feature, "_", view)) %>%
  # I average one more time over drug activity per line. This is necessary as D020 was imaged twice. In this particular case, I am creating the average activity score
  dplyr::group_by(sample, feature, view) %>% summarise(value = mean(value)) %>%
  arrange(feature, sample) %>%
  dplyr::select(sample, feature, value, view)

plot(mofa_morphology_agg$value, mofa_morphology_agg_new$value)
plot(mofa_morphology_sampled$value, mofa_morphology_sampled_new$value)
plot(mofa_morphology_sampled$value, mofa_morphology_agg$value)
plot(mofa_morphology_sampled_new$value, mofa_morphology_agg_new$value)

# new code without pre-mean
mofa_morphology_agg_new_premean <- readRDS(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_aggregate.Rds")) %>% 
  dplyr::filter(drug == "DMSO") %>% 
  mutate(rep = paste0("r", replicate)) %>% 
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  #filter(size_log > 7.5) %>%
  #group_by(line) %>% 
  #summarise_at(vars(contains("pc")), funs(mean)) %>% 
  #ungroup() %>% 
  dplyr::select(line, contains("pc")) %>%
  gather(pca, value, -line) %>% rename(feature = pca, sample = line) %>% mutate(view = "morphology_view") %>% 
  mutate(feature = paste0(feature, "_", view)) %>%
  # adding sort
  arrange(feature, sample)

x <- mofa_morphology_agg_new_premean %>% filter(sample != "D020_r1")  %>% filter(sample != "D020_r2")
y <- mofa_morphology_agg_new %>% filter(sample != "D020_r1")  %>% filter(sample != "D020_r2")
z <- mofa_morphology_agg %>% filter(sample != "D020_r1")  %>% filter(sample != "D020_r2")

data.frame(x$value,
       y$value, 
       z$value) %>% as.matrix() %>% cor()

# it is the group_by pre-mean calculation that messes it all up
