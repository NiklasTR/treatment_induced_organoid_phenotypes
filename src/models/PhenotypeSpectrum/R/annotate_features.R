library(rhdf5)
library(tidyverse)
library(here)

path = here("data/interim/FeatureAnalysis/line_differences/human/all_drugs/results/ReducedFeatures_all_drugs_human.h5")

print('extracting annotation')

h5ls(path)
df <- tibble(name = h5read(path, "feature_names_organoids")) %>%
  mutate(class = case_when(grepl(pattern = "\\.s\\.", x = name) ~ "shape",
                           grepl(pattern = "\\.b\\.", x = name) ~ "basic",
                           grepl(pattern = "\\.h\\.", x = name) ~ "haralick",
                           grepl(pattern = "\\.m\\.", x = name) ~ "moment")) %>%
  mutate(channel = case_when(grepl(pattern = "x\\.0\\.", x = name) ~ "mask",
                             grepl(pattern = "x\\.a\\.", x = name) ~ "actin",
                             grepl(pattern = "x\\.b\\.", x = name) ~ "cell_event",
                             grepl(pattern = "x\\.c\\.", x = name) ~ "dapi"))

# TODO build test to check wether the returned features equal the features listed in Comon_Features_human.txt in that directory

df %>% write_csv(here("data/interim/FeatureAnalysis/feature_annotation.csv"))

print('reading features')

df <- h5read(path, "features_organoids") %>% as.data.frame() %>% 
  magrittr::set_colnames(h5read(path, "feature_names_organoids")) %>%
  dplyr::select(x.a.b.mean, x.b.b.mean, x.c.b.mean)

df %>% write_rds(here("data/interim/FeatureAnalysis/feature_intensity.rds"))

print(head(df))

