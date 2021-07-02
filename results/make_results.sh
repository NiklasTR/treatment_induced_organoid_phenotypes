echo current working directory for R kernel:
Rscript -e "here::here()"
echo start rendering Rmarkdown vignettes and moving figures into subdirectories

# figure 1
echo figure 1
#Rscript -e "rmarkdown::render(here::here('notebooks/PhenotypeSpectrum/organoid_unsupervised_exploration.Rmd'))"
cp ../reports/panels/panel_size_dist.pdf figure_1/

# figure 2
cp ../reports/panels/panel_size_drug.pdf figure_2/