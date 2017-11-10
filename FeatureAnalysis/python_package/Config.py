import os

# Contains the hard-coded URLS so that these don't have to be passed between
# functions

BASEDIR = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/python_package"
# BASEDIR = "/home/sauerja/mnt/b210-projects/users/sauerja/PROMISE/python_package"
FEATUREDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
LAYOUTDIR = "/collab-ag-fischer/PROMISE/layouts/python_friendly"
SEGMENTATIONDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/segmentation"
FEATURETYPE = "organoids"
SIZETHRESHOLD = 1500
BLURRYORGANOIDDIR = os.path.join(BASEDIR, "blurry_organoids")
if not os.path.isdir(BLURRYORGANOIDDIR):
    os.makedirs(BLURRYORGANOIDDIR)
BLURRYWELLFN = os.path.join(
    BLURRYORGANOIDDIR, "blurry_wells_predicted.txt")
DEADORGANOIDCLASSIFIERDIR = os.path.join(BASEDIR, "dead_organoid_classifier")
DRUGTARGETCLASSIFIERDIR = os.path.join(BASEDIR, "drug_target_classifier")
