# packages
library(MOFA2)
library(tidyverse)
#library(here)

# input 
promise_long_filtered <- readRDS(here::here('data/processed/expression/promise_expr_filtered_tidy.rds'))

organoid_size_fit <- readRDS(here::here("data/processed/morphology/organoid_size.Rds")) %>% 
  filter(!line %in% c('D055T01', 'D020T02', 'D021T01')) %>% 
  #filter(!line %in% c('D055T01','D020T02')) %>% 
  mutate(line = as.character(line)) %>% 
  dplyr::select(line, size = x, rep = replicate) %>% 
  distinct() %>% arrange(line) %>%
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(rep = paste0("r", rep))

# define MOFA object and run MOFA
mofa_size <- organoid_size_fit %>% 
  group_by(line) %>% summarise(value = mean(size)) %>% 
  mutate(feature = "size",
         view = "size_view") %>% 
  drop_na() %>% 
  dplyr::select(sample = line, feature, view, value)

mofa_expression <- promise_long_filtered %>% 
  # renaming columns
  dplyr::select(sample = line,
                feature = probe,
                value = expr) %>% 
  mutate(view = "expression") %>% 
  # averaging feature value
  group_by(feature, sample, view) %>%
  summarise(value = mean(value)) %>%
  # renaming feature jic
  mutate(feature = paste0(feature, "_", view)) %>%
  dplyr::select(sample, feature, view, value) %>% 
  drop_na() # 

input_df = rbind(mofa_size,
                 mofa_expression) %>% 
  data.table::as.data.table()

MOFAobject <- create_mofa(input_df, verbose = TRUE)
print("created MOFA object")
plot_data_overview(MOFAobject) #%>% ggsave(here::here("reports/figures/mofa_object.pdf"))
system("mv Rplots.pdf reports/figures/mofa_object.pdf")

## setting option
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views = FALSE # default is TRUE
print(data_opts)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors = 5 # default is 15
print(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$verbose = TRUE
train_opts$maxiter = 1000 # 1000 is default
train_opts$stochastic <- FALSE # default FALSE
train_opts$save_interrupted <- TRUE
print(train_opts)

stochastic_opts <- get_default_stochastic_options(MOFAobject)
print(stochastic_opts)

# running model
outfile = file.path(here::here("models/mofa/model.hdf5"))
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  #stochastic_options = stochastic_opts, # comment out if not running stochastic inference
  training_options = train_opts
)
