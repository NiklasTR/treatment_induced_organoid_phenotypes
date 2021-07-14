# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
print(.libPaths())

# figure 1
rmarkdown::render(here::here('notebooks/organoid_unsupervised_exploration.Rmd'), params = list(remote = TRUE, data = "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_tidy.Rds", cache = FALSE))
#system("cp ../reports/panels/panel_size_dist.pdf figure_1/")

# figure 2
# cp ../reports/panels/panel_size_drug.pdf figure_2/
