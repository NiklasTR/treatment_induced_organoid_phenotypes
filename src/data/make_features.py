# exporting selected features from organoid feature table.
os.system('Rscript src/features/annotate_features.R')

# adding LDC classification results to UMAP projected (uncorrected and harmony corrected) data
os.system('Rscript src/models/PhenotypeSpectrum/R/LDC_organoids.R data/interim/FeatureAnalysis/organoid_viability/human/results data/interim/PhenotypeSpectrum/harmony_umap_absolute_all_drugs.Rds')
os.system('Rscript src/models/PhenotypeSpectrum/R/LDC_organoids.R data/interim/FeatureAnalysis/organoid_viability/human/results data/interim/PhenotypeSpectrum/hdf5_umap_absolute_all_drugs.Rds')

# exporting tidied and small subsets of the data for interactive handling.
os.system('Rscript src/features/annotate_features.R')

os.system('Rscript src/data/make_expression.R')
