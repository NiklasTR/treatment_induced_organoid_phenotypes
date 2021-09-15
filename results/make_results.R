# defining lib path
#.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
print(.libPaths())

# figure 1
rmarkdown::render(here::here('notebooks/imaging/organoid_unsupervised_exploration.Rmd'), 
    params = list(
        remote = FALSE, 
        data = "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds", 
        sample = "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_tidy_Paclitaxel.Rds",
        cache = TRUE))

# figure 2
rmarkdown::render(here::here('notebooks/imaging/embedding_inspection.Rmd'), 
    params = list(
        data = "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds",
        data_harmony = "data/processed/PhenotypeSpectrum/harmony_umap_absolute_all_drugs_sampled.Rds",
        remote = FALSE,
        cache = TRUE))
rmarkdown::render(here::here('notebooks/drug_analysis/2.0-js-OrganoidViability.Rmd'))

# figure 3
rmarkdown::render(here::here('notebooks/drug_analysis/3.0-js-OrganoidPhenotypes.Rmd'))

# figure 4 and 5
rmarkdown::render(here::here('notebooks/mofa/mofa_exploration.Rmd'))
rmarkdown::render(here::here('notebooks/mofa/mofa_drug_effect.Rmd'))