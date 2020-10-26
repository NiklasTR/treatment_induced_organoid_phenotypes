import os
######################################
# ONLY runs on system with IBM Cluster, e.g. ssh rindtorf@odcf-lsf01.dkfz.de
# cd /dkfz/groups/shared/OE0049/B110-Isilon2/promise
print('make sure you are running make_dataset on an instance that supports IBM bsub')

# running incremental PCA on data
os.system('bsub -R "rusage[mem=150GB]" -q long ./src/features/FeatureAnalysis/run_pca.bsub \
	-o src/features/FeatureAnalysis/out.bsub \
	-e src/features/FeatureAnalysis/error.bsub')
# Aggregated features 6.6M x 486 and incrementalPCA projected features are stored in:  
# data/interim/FeatureAnalysis/line_differences/human/all_drugs/results

# 
os.system('bsub -R "rusage[mem=100GB]" -q long ./src/models/PhenotypeSpectrum/run_umap_dmso.bsub \
	-o src/models/PhenotypeSpectrum/out.bsub \
	-e src/models/PhenotypeSpectrum/error.bsub')
	
os.system('bsub -R "rusage[mem=200GB]" -q verylong ./src/models/PhenotypeSpectrum/run_umap_all_drugs.bsub \
	-o src/models/PhenotypeSpectrum/out.bsub \
	-e src/models/PhenotypeSpectrum/error.bsub')


# exporting selected features from organoid feature table.
os.system('Rscript src/features/annotate_features.R')




# adding LDC classification results to UMAP projected (uncorrected and harmony corrected) data
os.system('Rscript src/models/PhenotypeSpectrum/R/LDC_organoids.R data/interim/FeatureAnalysis/organoid_viability/human/results data/interim/PhenotypeSpectrum/harmony_umap_absolute_all_drugs.Rds')
os.system('Rscript src/models/PhenotypeSpectrum/R/LDC_organoids.R data/interim/FeatureAnalysis/organoid_viability/human/results data/interim/PhenotypeSpectrum/hdf5_umap_absolute_all_drugs.Rds')
