library(here)
library(tidyverse)
library(MOFA2)
library(data.table)


mofa_df <- read_rds(here("data/processed/PhenotypeSpectrum/mofa_df.Rds"))

# brute force reduction
set.seed(123)
mofa_df_filtered <- mofa_df %>% 
  mutate(sample = paste0(sample, "_", group)) %>% 
  mutate(feature = paste0(feature, "_", view)) %>% 
  # log transforming size
  mutate(value = ifelse(feature == "x.0.s.area_expected", log(value), value)) %>%
  filter(view != "mutation") %>%
  # scaling features manually
  group_by(feature, view) %>% 
  mutate(value = if_else(view == "morphology", scale(value), scale(value, center = FALSE))) %>%
  mutate(value = scale(value)) %>%
  ungroup() 

mofa_df_filtered %>% write_rds(here("data/processed/PhenotypeSpectrum/mofa_df_filtered.Rds"))

input_df <- mofa_df_filtered %>%
  filter(grepl(sample, pattern = "DMSO")) %>%
  # mutate(feature = paste0(feature, "_", view)) %>% 
  dplyr::select(-group) %>% 
  as.data.table() %>% 
  drop_na()
  
print("wrangled data")
MOFAobject <- create_mofa(input_df, verbose = TRUE)
print("created MOFA object")
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views = TRUE # default is TRUE
print(data_opts)

model_opts <- get_default_model_options(MOFAobject)
print(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$verbose = TRUE
train_opts$maxiter = 1000 # 1000 is default
train_opts$stochastic <- FALSE
train_opts$save_interrupted <- TRUE
print(train_opts)

stochastic_opts <- get_default_stochastic_options(MOFAobject)
print(stochastic_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  #stochastic_options = stochastic_opts,
  training_options = train_opts
)


outfile = file.path(here("data/processed/PhenotypeSpectrum/mofa_model_dmso.hdf5"))
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

