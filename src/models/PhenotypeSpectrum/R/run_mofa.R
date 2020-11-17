library(here)
library(tidyverse)
library(MOFA2)
library(data.table)


mofa_df <- read_rds(here("data/processed/PhenotypeSpectrum/mofa_df.Rds"))

# brute force reduction
set.seed(123)
input_df <- mofa_df %>% mutate(sample = paste0(sample, "_", group)) %>% 
  dplyr::select(-group) %>% 
  #sample_frac(1) %>% 
  filter(view != "mutation") %>%
  as.data.table()
 
  
print("wrangled data")
MOFAobject <- create_mofa(input_df, verbose = TRUE)
print("created MOFA object")
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views = TRUE
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$verbose = TRUE
head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)


outfile = file.path(here("data/processed/PhenotypeSpectrum/mofa_model.hdf5"))
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

