echo current working directory for R kernel:
Rscript -e "here::here()"
Rscript -e ".libPaths()"
Rscript -e ".libPaths(here::here("promise/x86_64-pc-linux-gnu-library/3.6"))"
echo start rendering Rmarkdown vignettes and moving figures into subdirectories

# figure 1
echo figure 1
Rscript -e "rmarkdown::render(here::here('notebooks/organoid_unsupervised_exploration.Rmd'), params = list(remote = TRUE))"
#cp ../reports/panels/panel_size_dist.pdf figure_1/

# figure 2
# cp ../reports/panels/panel_size_drug.pdf figure_2/
