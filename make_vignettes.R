library(tidyverse)

tibble(vignettes = list.files(here::here("vignettes"), pattern = ".Rmd", full.names = TRUE) %>% sort()) %>%
  mutate(scripts = vignettes %>% str_replace(pattern = ".Rmd", replacement = ".R")) %>% 
  mutate(purl = map2(vignettes, scripts, ~ .x %>% knitr::purl(output = .y))) %>% 
  mutate(render = walk(vignettes, ~ .x %>% rmarkdown::render())) # %>%
 # mutate(source = walk(scripts, ~ .x %>% source()))
  