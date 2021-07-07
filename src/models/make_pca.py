import os
######################################
# ONLY runs on system with IBM Cluster, e.g. ssh rindtorf@odcf-lsf01.dkfz.de
# cd /dkfz/groups/shared/OE0049/B110-Isilon2/promise
print('make sure you are running make_dataset on an instance that supports IBM bsub')

# 1. running incremental PCA on data and exporting an .Rds object
# ====================================
os.system('bsub -R "rusage[mem=150GB]" -q long ./src/models/FeatureAnalysis/run_pca.bsub \
	-o src/models/FeatureAnalysis/run_pca_out.bsub \
	-e src/models/FeatureAnalysis/run_pca_error.bsub')

