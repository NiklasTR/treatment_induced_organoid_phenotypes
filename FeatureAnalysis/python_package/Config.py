import os

# Contains the hard-coded URLS so that these don't have to be passed between
# functions

BASEDIR = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/python_package"
# BASEDIR = "/home/sauerja/mnt/b210-projects/users/sauerja/PROMISE/preprocessing"
FEATUREDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
LAYOUTDIR = "/collab-ag-fischer/PROMISE/layouts/python_friendly"
SEGMENTATIONDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/segmentation"
BLURRYWELLFN = os.path.join(BASEDIR, "blurry_wells_predicted.txt")
FEATURETYPE = "organoids"
SIZETHRESHOLD = 2500
BLURRYORGANOIDCLF = os.path.join(
    BASEDIR, "blurry_organoid_classifier_%s.pkl" % FEATURETYPE)
