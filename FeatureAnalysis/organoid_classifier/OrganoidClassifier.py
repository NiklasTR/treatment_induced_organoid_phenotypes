# This script classifies individual organoids from negative and positive
# controls

from __future__ import division
import os
import random
import h5py
import numpy as np
import sklearn.model_selection
import sklearn.ensemble
import re
import pickle
import statsmodels.robust
import pandas as pd
import scipy.stats.mstats
import sys

BASEDIR = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/organoid_classifier"
FEATUREDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
LAYOUTDIR = "/collab-ag-fischer/PROMISE/layouts/python_friendly"
# SEGMENTATIONDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/segmentation"
# BLURRYWELLFN = os.path.join(BASEDIR, "blurry_wells_predicted.txt")
FEATURETYPE = "organoids"
SIZETHRESHOLD = 2500


def get_wells_for_treatment(plate, treatment):
    """
    Get the well IDs corresponding to a given treatment
    :return:
    """

    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(LAYOUTDIR, "%s.xlsx" % layout_id))
    wells = layout.loc[
        layout["Product.Name"] == "DMSO",
        "Well_ID_384"].values
    return wells


def
