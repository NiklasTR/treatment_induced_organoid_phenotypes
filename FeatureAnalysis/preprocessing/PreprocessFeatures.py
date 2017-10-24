# This is the feature preprocessing script
# The general steps are:
#   - Train and apply a classifier to learn the difference between blurry
#     and in-focus organoids.
#   - Calculate an average for each well
#     - Filter out blurry organoids
#     - Separate into organoids and shrapnel (currently with threshold)
#     - Calc trimmed mean and sd
#   - Combine well features into plate feature

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


BASEDIR = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/preprocessing"
FEATUREDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
LAYOUTDIR = "/collab-ag-fischer/PROMISE/layouts/python_friendly"
SEGMENTATIONDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/segmentation"
BLURRYWELLFN = os.path.join(BASEDIR, "blurry_wells_predicted.txt")
FEATURETYPE = "organoids"
SIZETHRESHOLD = 2500


def trim_func(a, func, percent, axis):
    """
    Apply a function to a trimmed matrix along an axis. The function must be
    of a type that takes the 'axis' keyword, e.g. np.mean, np.std, np.median,
    etc.
    :param func:
    :param a:
    :param percent:
    :param axis:
    :return:
    """
    asort = np.sort(a, axis=axis)
    num_cases = asort.shape[1]
    num_to_trim = int(np.ceil(percent*num_cases))
    atrim = np.delete(
        asort, slice(num_cases - num_to_trim, num_cases),
        axis=axis)
    atrim = np.delete(
        atrim, slice(0, num_to_trim),
        axis=axis)
    if atrim.shape[1] == 0:
        return np.zeros(shape=a.shape[0])
    else:
        return func(atrim, axis=axis)


def learn_blurry_organoids():
    """
    Trains a classifier to tell the difference between blurry and in-focus
    organoids based on the features. This function currently explicitly
    doesn't use the blurry wells detected via inception features. This may
    be unnecessary but why risk it.

    The steps are:
    1. Initial annotation with a high-precision/low-recall segmentation based
       on maximal DAPI intensity.
    2. Train random forest classifier on all unrelated features to prevent
       overfitting

    :return:
        - A random forest classifier object
        - The names of features used to train the RF
        - The expected accuracy on a validation set
    """

    keymap = {
        "organoids": ("features_organoids", "feature_names_organoids"),
        "clumps": ("features_clumps", "feature_names_clumps")}

    if FEATURETYPE not in keymap.keys():
        raise KeyError("'feature_type' must be one of '%s'"
                       % str(keymap.keys()))

    hdf5_keys = keymap[FEATURETYPE]

    with open(BLURRYWELLFN, "r") as f:
        blurry_wells = [s.strip() for s in f.readlines()]
    all_plates = [
        s for s in os.listdir(FEATUREDIR) if
        s.startswith("M001") or s.startswith("D0")]
    features = []
    feature_names = []
    for plate in all_plates:
        wells = [s for s in os.listdir(
            os.path.join(FEATUREDIR, plate, "wells"))]
        wells = [s for s in wells if s[0:19] not in blurry_wells]
        wells = random.sample(wells, 15)
        for well in wells:
            feature_fn = os.path.join(FEATUREDIR, plate, "wells", well)
            try:
                with h5py.File(feature_fn, "r") as h5handle:
                    features.append(h5handle[hdf5_keys[0]][()])
                    feature_names.append(h5handle[hdf5_keys[1]][()])
            except KeyError:
                pass
    features = np.concatenate(features, axis=1)
    features = features.transpose()
    feature_names = feature_names[0]

    features = np.delete(
        features, np.where(feature_names == "FIELD")[0], axis=1)
    feature_names = np.delete(
        feature_names, np.where(feature_names == "FIELD")[0])

    cy3_ind = np.where(feature_names == "x.a.b.q099")[0][0]
    max_cy3_intensity = features[:, cy3_ind]
    max_cy3_thresh = np.percentile(max_cy3_intensity, 90)
    min_cy3_thresh = np.percentile(max_cy3_intensity, 35)
    dapi_ind = np.where(feature_names == "x.c.b.q099")[0][0]
    max_dapi_intensity = features[:, dapi_ind]
    max_dapi_thresh = np.percentile(max_dapi_intensity, 90)
    min_dapi_thresh = np.percentile(max_dapi_intensity, 35)
    pos_training = features[np.logical_and(
        max_dapi_intensity >= max_dapi_thresh,
        max_cy3_thresh >= max_cy3_thresh), :]
    neg_training = features[np.logical_and(
        max_dapi_intensity <= min_dapi_thresh,
        max_cy3_intensity <= min_cy3_thresh), :]
    min_size = min(pos_training.shape[0], neg_training.shape[0])
    pos_rand_ind = np.random.choice(
        a=range(pos_training.shape[0]), size=min_size)
    pos_training = pos_training[pos_rand_ind, :]
    neg_rand_ind = np.random.choice(
        a=range(neg_training.shape[0]), size=min_size)
    neg_training = neg_training[neg_rand_ind, :]

    X = np.concatenate((neg_training, pos_training), axis=0)
    Y = np.repeat((0, 1), min_size)

    # Remove the feature used to annotate, and similar features, to prevent
    # overfitting
    remove_features_names = [
        s for s in feature_names if
        re.match("x\..*\.b", s)
        is not None]
    altered_feature_names = [
        s for s in feature_names if s not in remove_features_names]
    index_mask = np.ma.array(range(len(feature_names)), mask=False)
    for feat_name in remove_features_names:
        feat_index = np.where(feature_names == feat_name)[0][0]
        index_mask[feat_index] = np.ma.masked

    X = X[:, index_mask.compressed()]
    X_train, X_val, Y_train, Y_val = sklearn.model_selection.train_test_split(
        X, Y, test_size=0.25)

    rfclf = sklearn.ensemble.RandomForestClassifier(n_estimators=100)
    rfclf.fit(X_train, Y_train)
    val_acc = rfclf.score(X_val, Y_val)

    return {
        "clf": rfclf, "feature_names": altered_feature_names,
        "accuracy": val_acc}


def calc_well_average(
        features, feature_names, blurry_organoid_clf,
        summary_func_middle, summary_func_var, kwargs_middle,
        kwargs_var):
    """
    Calculate the summaries over features.

    'blurry_well_clf' should be the dictionary object returned by the
    function 'learn_blurry_organoids()'

    :param features: A 2D numpy matrix with the shape (features, samples)
    :param feature_names:
    :param blurry_organoid_clf: A dictionary.
    :param summary_func_middle:
    :param summary_func_var:
    :param kwargs_middle:
    :param kwargs_var:

    :return:
        - A feature vector for the well
        - A list of feature names
    """

    # Remove the feature "FIELD"
    features = np.delete(
        features, np.where(feature_names == "FIELD")[0], axis=0)
    feature_names = np.delete(
        feature_names, np.where(feature_names == "FIELD")[0])

    features_clf = np.array(features).transpose()
    features_clf = features_clf[
        ..., [f in blurry_organoid_clf["feature_names"]
              for f in feature_names]]

    is_focused = blurry_organoid_clf["clf"].predict(features_clf)
    features = features[..., [s == 1 for s in is_focused]]

    # Separate into organoids and shrapnel
    f_size = features[np.where(feature_names == "x.0.s.area")[0][0]]
    features_shrapnel = features[:, f_size < SIZETHRESHOLD]
    features_organoids = features[:, f_size >= SIZETHRESHOLD]

    # Calculate feature summaries (median)
    if features_shrapnel.shape[1] == 0:
        features_shrapnel_m = np.zeros(
            shape=features_shrapnel.shape[0])
        features_shrapnel_v = np.zeros(
            shape=features_shrapnel.shape[0])
    else:
        features_shrapnel_m = summary_func_middle(
            features_shrapnel, axis=1, **kwargs_middle)
        features_shrapnel_v = summary_func_var(
            features_shrapnel, axis=1, **kwargs_var)
    features_shrapnel_summary = np.concatenate(
        (features_shrapnel_m, features_shrapnel_v,
         (features_shrapnel.shape[1],)))
    
    if features_organoids.shape[1] == 0:
        features_organoids_m = np.zeros(
            shape=features_organoids.shape[0])
        features_organoids_v = np.zeros(
            shape=features_organoids.shape[0])
    else:
        features_organoids_m = summary_func_middle(
            features_organoids, axis=1, **kwargs_middle)
        features_organoids_v = summary_func_var(
            features_organoids, axis=1, **kwargs_var)
    features_organoids_summary = np.concatenate(
        (features_organoids_m, features_organoids_v,
         (features_organoids.shape[1],)))

    features_summary = np.concatenate(
        (features_organoids_summary, features_shrapnel_summary))

    full_feature_names = np.array(
        ["organoids_%s_expected" % fn for fn in feature_names] +
        ["organoids_%s_variation" % fn for fn in feature_names] +
        ["organoids_num.of.objects"] +
        ["shrapnel_%s_expected" % fn for fn in feature_names] +
        ["shrapnel_%s_variation" % fn for fn in feature_names] +
        ["shrapnel_num.of.objects"])

    return features_summary, full_feature_names


def transform_features(features, transform_type, **kwargs):
    """
    Performs a transformation on the features
    :param features: A 2D numpy array with the shape (features, samples)
    :param transform_type: The type of transformation to perform
    :param kwargs: Keyword arguments to the transformation type
    :return:
    """

    if transform_type == "boxcox":
        raise Warning("Box-Cox transformation not implemented yet")
        # # Values must first be pushed into the positive
        # fmin = np.min(f, axis=1)
        # tmp2 = np.apply_along_axis(func1d=scipy.stats.boxcox, axis=axis, arr=f)
        # norm_features = scipy.stats.boxcox()
    elif transform_type == "glog":
        if "c" in kwargs.keys():
            c = kwargs["c"]
        else:
            c = 0.05
        return np.log((features + np.sqrt(features ** 2 + c ** 2)) / 2.0)
    else:
        raise ValueError("Unsupported transform_type ('%s')" % transform_type)


def calc_plate_median(plate, blurry_organoid_clf):
    """
    Calculate and return the plate summary.

    It calculates:
    - median / mad
    - trimmed mean / trimmed sd @ 5% trimming
    - trimmed mean / trimmed sd @ 10% trimming
    - trimmed mean / trimmed sd @ 15% trimming

    :param plate:
    :param blurry_organoid_clf:
    :return:
    """
    wells = sorted([
        well for well in
        os.listdir(os.path.join(FEATUREDIR, plate, "wells"))
        if well.startswith(plate) if well.endswith(".h5")])
    median_features = []
    trimmed_features_05 = []
    trimmed_features_10 = []
    trimmed_features_15 = []
    plate_feature_names = []
    for well in wells:
        well_fn = os.path.join(FEATUREDIR, plate, "wells", well)
        try:
            with h5py.File(well_fn, "r") as h5handle:
                features = h5handle["features_%s" % FEATURETYPE][()]
                feature_names = h5handle["feature_names_%s" % FEATURETYPE][()]
            median = calc_well_average(
                features=features, feature_names=feature_names,
                blurry_organoid_clf=blurry_organoid_clf,
                summary_func_middle=np.nanmedian,
                summary_func_var=statsmodels.robust.mad,
                kwargs_middle=dict(), kwargs_var={"center": np.nanmedian})
            trimmed_05 = calc_well_average(
                features=features, feature_names=feature_names,
                blurry_organoid_clf=blurry_organoid_clf,
                summary_func_middle=trim_func,
                summary_func_var=trim_func,
                kwargs_middle={"func": np.mean, "percent": 0.05},
                kwargs_var={"func": np.std, "percent": 0.05})
            trimmed_10 = calc_well_average(
                features=features, feature_names=feature_names,
                blurry_organoid_clf=blurry_organoid_clf,
                summary_func_middle=trim_func,
                summary_func_var=trim_func,
                kwargs_middle={"func": np.mean, "percent": 0.10},
                kwargs_var={"func": np.std, "percent": 0.10})
            trimmed_15 = calc_well_average(
                features=features, feature_names=feature_names,
                blurry_organoid_clf=blurry_organoid_clf,
                summary_func_middle=trim_func,
                summary_func_var=trim_func,
                kwargs_middle={"func": np.mean, "percent": 0.15},
                kwargs_var={"func": np.std, "percent": 0.15})
            if not np.all(median[1] == trimmed_05[1]):
                raise Exception("Feature names not identical")
            if not np.all(median[1] == trimmed_10[1]):
                raise Exception("Feature names not identical")
            if not np.all(median[1] == trimmed_15[1]):
                raise Exception("Feature names not identical")
            median_features.append(median[0])
            trimmed_features_05.append(trimmed_05[0])
            trimmed_features_10.append(trimmed_10[0])
            trimmed_features_15.append(trimmed_15[0])
            plate_feature_names.append(median[1])
        except Exception:
            median_features.append(np.repeat(np.nan, 6286))
            trimmed_features_05.append(np.repeat(np.nan, 6286))
            trimmed_features_10.append(np.repeat(np.nan, 6286))
            trimmed_features_15.append(np.repeat(np.nan, 6286))

    fname_iter = iter(plate_feature_names)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    plate_feature_names = plate_feature_names[0]

    median_features = np.stack(median_features, axis=1)
    trimmed_features_05 = np.stack(trimmed_features_05, axis=1)
    trimmed_features_10 = np.stack(trimmed_features_10, axis=1)
    trimmed_features_15 = np.stack(trimmed_features_15, axis=1)
    well_names = np.array([well[0:19] for well in wells])
    return {
        "features_median": median_features,
        "features_trimmed_mean_05": trimmed_features_05,
        "features_trimmed_mean_10": trimmed_features_10,
        "features_trimmed_mean_15": trimmed_features_15,
        "feature_names": plate_feature_names,
        "well_names": well_names}


def calc_dmso_avg_of_plate(plate, blurry_organoid_clf):
    """
    Calculate the DMSO average (median) of a plate
    :param plate:
    :param blurry_organoid_clf:
    :return:
    """
    out_fn = os.path.join(
        FEATUREDIR, plate,
        "%s_dmso_summaries_%s.h5" % (plate, FEATURETYPE))
    if os.path.isfile(out_fn):
        return False

    # Load the library of the plate and extract the DMSO wells names
    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(LAYOUTDIR, "%s.xlsx" % layout_id))
    dmso_wells = layout.loc[
        layout["Product.Name"] == "DMSO",
        "Well_ID_384"].values

    # Load the individual organoid features of the wells
    dmso_features = []
    dmso_feature_names = []
    num_organoids = []
    num_shrapnel = []
    for dmso_well in dmso_wells:
        well_fn = os.path.join(
            FEATUREDIR, plate, "wells",
            "%s_%s_%s_features.h5" %
            (plate, dmso_well[0], dmso_well[1:3]))
        try:
            with h5py.File(well_fn, "r") as h5handle:
                well_features = h5handle["features_%s" % FEATURETYPE][()]
                well_feature_names = h5handle["feature_names_%s" % FEATURETYPE][()]
                well_f_size = well_features[np.where(
                    well_feature_names == "x.0.s.area")[0][0]]
                num_organoids.append(np.sum(well_f_size >= SIZETHRESHOLD))
                num_shrapnel.append(np.sum(well_f_size < SIZETHRESHOLD))
                dmso_features.append(well_features)
                dmso_feature_names.append(well_feature_names)
        except KeyError:
            pass
    dmso_features = np.concatenate(dmso_features, axis=1)
    fname_iter = iter(dmso_feature_names)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    dmso_feature_names = dmso_feature_names[0]

    # Calculate the summaries over the features
    median = calc_well_average(
        features=dmso_features, feature_names=dmso_feature_names,
        blurry_organoid_clf=blurry_organoid_clf,
        summary_func_middle=np.nanmedian,
        summary_func_var=statsmodels.robust.mad,
        kwargs_middle=dict(), kwargs_var={"center": np.nanmedian})
    trimmed_05 = calc_well_average(
        features=dmso_features, feature_names=dmso_feature_names,
        blurry_organoid_clf=blurry_organoid_clf,
        summary_func_middle=trim_func,
        summary_func_var=trim_func,
        kwargs_middle={"func": np.mean, "percent": 0.05},
        kwargs_var={"func": np.std, "percent": 0.05})
    trimmed_10 = calc_well_average(
        features=dmso_features, feature_names=dmso_feature_names,
        blurry_organoid_clf=blurry_organoid_clf,
        summary_func_middle=trim_func,
        summary_func_var=trim_func,
        kwargs_middle={"func": np.mean, "percent": 0.10},
        kwargs_var={"func": np.std, "percent": 0.10})
    trimmed_15 = calc_well_average(
        features=dmso_features, feature_names=dmso_feature_names,
        blurry_organoid_clf=blurry_organoid_clf,
        summary_func_middle=trim_func,
        summary_func_var=trim_func,
        kwargs_middle={"func": np.mean, "percent": 0.15},
        kwargs_var={"func": np.std, "percent": 0.15})
    if not np.all(median[1] == trimmed_05[1]):
        raise Exception("Feature names not identical")
    if not np.all(median[1] == trimmed_10[1]):
        raise Exception("Feature names not identical")
    if not np.all(median[1] == trimmed_15[1]):
        raise Exception("Feature names not identical")

    # I'm "misusing" the function calc_well_average. I need to separate the
    # expected and deviation values for each entry as well as calculate the
    # proper expected value and variation for the number of objects
    def separate_featset(fs, summary_type):
        summary_dict = {
            "median": (np.nanmedian, statsmodels.robust.mad,
                       dict(), {"center": np.nanmedian}),
            "trimmed_05": (trim_func, trim_func,
                           {"func": np.mean, "percent": 0.05},
                           {"func": np.std, "percent": 0.05}),
            "trimmed_10": (trim_func, trim_func,
                           {"func": np.mean, "percent": 0.10},
                           {"func": np.std, "percent": 0.10}),
            "trimmed_15": (trim_func, trim_func,
                           {"func": np.mean, "percent": 0.15},
                           {"func": np.std, "percent": 0.15})}
        summary_func = summary_dict[summary_type]
        num_organoids_np = np.expand_dims(num_organoids, axis=0)
        num_shrapnel_np = np.expand_dims(num_shrapnel, axis=0)

        fs_org_exp_index = np.array([
            re.match(pattern="^organoids_.*_expected$", string=featname)
            is not None for featname in fs[1]])
        org_exp_feat = fs[0][fs_org_exp_index]
        org_exp_feat_name = fs[1][fs_org_exp_index]
        num_organoids_expected = summary_func[0](
            num_organoids_np, axis=1, **summary_func[2])
        org_exp_feat = np.concatenate((org_exp_feat, num_organoids_expected))
        org_exp_feat_name = np.append(
            org_exp_feat_name, "organoids_num.of.objects_expected")

        fs_org_dev_index = np.array([
            re.match(pattern="^organoids_.*_variation$", string=featname)
            is not None for featname in fs[1]])
        org_dev_feat = fs[0][fs_org_dev_index]
        org_dev_feat_name = fs[1][fs_org_dev_index]
        num_organoids_dev = summary_func[1](
            num_organoids_np, axis=1, **summary_func[3])
        org_dev_feat = np.concatenate((org_dev_feat, num_organoids_dev))
        org_dev_feat_name = np.append(
            org_dev_feat_name, "organoids_num.of.objects_variation")

        fs_shr_exp_index = np.array([
            re.match(pattern="^shrapnel_.*_expected$", string=featname)
            is not None for featname in fs[1]])
        shr_exp_feat = fs[0][fs_shr_exp_index]
        shr_exp_feat_name = fs[1][fs_shr_exp_index]
        num_shrapnel_expected = summary_func[0](
            num_shrapnel_np, axis=1, **summary_func[2])
        shr_exp_feat = np.concatenate((shr_exp_feat, num_shrapnel_expected))
        shr_exp_feat_name = np.append(
            shr_exp_feat_name, "shrapnel_num.of.objects_expected")

        fs_shr_dev_index = np.array([
            re.match(pattern="^shrapnel_.*_variation$", string=featname)
            is not None for featname in fs[1]])
        shr_dev_feat = fs[0][fs_shr_dev_index]
        shr_dev_feat_name = fs[1][fs_shr_dev_index]
        num_shrapnel_dev = summary_func[1](
            num_shrapnel_np, axis=1, **summary_func[3])
        shr_dev_feat = np.concatenate((shr_dev_feat, num_shrapnel_dev))
        shr_dev_feat_name = np.append(
            shr_dev_feat_name, "shrapnel_num.of.objects_variation")

        return {"organoids_expected": (org_exp_feat, org_exp_feat_name),
                "organoids_variation": (org_dev_feat, org_dev_feat_name),
                "shrapnel_expected": (shr_exp_feat, shr_exp_feat_name),
                "shrapnel_variation": (shr_dev_feat, shr_dev_feat_name)}

    median = separate_featset(median, "median")
    trimmed_05 = separate_featset(trimmed_05, "trimmed_05")
    trimmed_10 = separate_featset(trimmed_10, "trimmed_10")
    trimmed_15 = separate_featset(trimmed_15, "trimmed_15")

    with h5py.File(out_fn, "w") as h5handle:
        h5handle.create_dataset(
            name="dmso_median_organoids",
            data=median["organoids_expected"][0])
        h5handle.create_dataset(
            name="dmso_median_organoids_names",
            data=median["organoids_expected"][1])
        h5handle.create_dataset(
            name="dmso_median_shrapnel",
            data=median["shrapnel_expected"][0])
        h5handle.create_dataset(
            name="dmso_median_shrapnel_names",
            data=median["shrapnel_expected"][1])
        h5handle.create_dataset(
            name="dmso_mad_organoids",
            data=median["organoids_variation"][0])
        h5handle.create_dataset(
            name="dmso_mad_organoids_names",
            data=median["organoids_variation"][1])
        h5handle.create_dataset(
            name="dmso_mad_shrapnel",
            data=median["shrapnel_variation"][0])
        h5handle.create_dataset(
            name="dmso_mad_shrapnel_names",
            data=median["shrapnel_variation"][1])

        h5handle.create_dataset(
            name="dmso_tm05_mean_organoids",
            data=trimmed_05["organoids_expected"][0])
        h5handle.create_dataset(
            name="dmso_tm05_mean_organoids_names",
            data=trimmed_05["organoids_expected"][1])
        h5handle.create_dataset(
            name="dmso_tm05_mean_shrapnel",
            data=trimmed_05["shrapnel_expected"][0])
        h5handle.create_dataset(
            name="dmso_tm05_mean_shrapnel_names",
            data=trimmed_05["shrapnel_expected"][1])
        h5handle.create_dataset(
            name="dmso_tm05_std_organoids",
            data=trimmed_05["organoids_variation"][0])
        h5handle.create_dataset(
            name="dmso_tm05_std_organoids_names",
            data=trimmed_05["organoids_variation"][1])
        h5handle.create_dataset(
            name="dmso_tm05_std_shrapnel",
            data=trimmed_05["shrapnel_variation"][0])
        h5handle.create_dataset(
            name="dmso_tm05_std_shrapnel_names",
            data=trimmed_05["shrapnel_variation"][1])

        h5handle.create_dataset(
            name="dmso_tm10_mean_organoids",
            data=trimmed_10["organoids_expected"][0])
        h5handle.create_dataset(
            name="dmso_tm10_mean_organoids_names",
            data=trimmed_10["organoids_expected"][1])
        h5handle.create_dataset(
            name="dmso_tm10_mean_shrapnel",
            data=trimmed_10["shrapnel_expected"][0])
        h5handle.create_dataset(
            name="dmso_tm10_mean_shrapnel_names",
            data=trimmed_10["shrapnel_expected"][1])
        h5handle.create_dataset(
            name="dmso_tm10_std_organoids",
            data=trimmed_10["organoids_variation"][0])
        h5handle.create_dataset(
            name="dmso_tm10_std_organoids_names",
            data=trimmed_10["organoids_variation"][1])
        h5handle.create_dataset(
            name="dmso_tm10_std_shrapnel",
            data=trimmed_10["shrapnel_variation"][0])
        h5handle.create_dataset(
            name="dmso_tm10_std_shrapnel_names",
            data=trimmed_10["shrapnel_variation"][1])

        h5handle.create_dataset(
            name="dmso_tm15_mean_organoids",
            data=trimmed_15["organoids_expected"][0])
        h5handle.create_dataset(
            name="dmso_tm15_mean_organoids_names",
            data=trimmed_15["organoids_expected"][1])
        h5handle.create_dataset(
            name="dmso_tm15_mean_shrapnel",
            data=trimmed_15["shrapnel_expected"][0])
        h5handle.create_dataset(
            name="dmso_tm15_mean_shrapnel_names",
            data=trimmed_15["shrapnel_expected"][1])
        h5handle.create_dataset(
            name="dmso_tm15_std_organoids",
            data=trimmed_15["organoids_variation"][0])
        h5handle.create_dataset(
            name="dmso_tm15_std_organoids_names",
            data=trimmed_15["organoids_variation"][1])
        h5handle.create_dataset(
            name="dmso_tm15_std_shrapnel",
            data=trimmed_15["shrapnel_variation"][0])
        h5handle.create_dataset(
            name="dmso_tm15_std_shrapnel_names",
            data=trimmed_15["shrapnel_variation"][1])

    return True


def normalize_organoids_in_well(well_id, blurry_organoid_clf):
    """
    Normalize all the individual organoid features in a well.
    This includes the same steps as for the well-averaged features
    but without the initial averaging step:

    - Filter out blurry organoids
    - Split into organoids and shrapnel
    - glog transform features
    - Normalize by DMSO values (f = (f - dmso_expected) / dmso_var)

    :param well_id:
    :param blurry_organoid_clf:
    :return:
    """
    plate, row, col = well_id.split("_")
    in_fn = os.path.join(
        FEATUREDIR, plate, "wells",
        "%s_%s_%s_features.h5" % (plate, row, col))
    out_fn = os.path.join(
        FEATUREDIR, plate, "wells_normalized"
        "%s_%s_%s_features_normalized.h5" % (plate, row, col))
    if not os.path.isdir(os.path.join(FEATUREDIR, plate, "wells_normalized")):
        os.makedirs(os.path.join(FEATUREDIR, plate, "wells_normalized"))
    if os.path.isfile(out_fn):
        return False

    # Load features for organoids and clumps
    features_dict = dict()
    feature_names_dict = dict()
    try:
        with h5py.File(in_fn, "r") as h5handle:
            features_dict["organoids"] = h5handle["features_organoids"][()]
            feature_names_dict["organoids"] = h5handle[
                "feature_names_organoids"][()]
            features_dict["clumps"] = h5handle["features_clumps"][()]
            feature_names_dict["clumps"] = h5handle[
                "feature_names_clumps"][()]
    except KeyError:
        return False

    # Loop the entire workflow over both keywords
    for featuretype in features_dict.keys():
        features = features_dict[featuretype]
        feature_names = feature_names_dict[featuretype]

        # Remove FIELD info
        features = np.delete(
            features, np.where(feature_names == "FIELD")[0], axis=0)
        feature_names = np.delete(
            feature_names, np.where(feature_names == "FIELD")[0])

        # Filter out blurry organoids
        features_clf = np.array(features).transpose()
        features_clf = features_clf[
            ..., [f in blurry_organoid_clf["feature_names"]
                  for f in feature_names]]

        is_focused = blurry_organoid_clf["clf"].predict(features_clf)
        features = features[..., [s == 1 for s in is_focused]]

        # Separate into organoids and shrapnel
        f_size = features[np.where(
            feature_names == "x.0.s.area")[0][0]]
        features_shrapnel = features[:, f_size < SIZETHRESHOLD]
        features_organoids = features[:, f_size >= SIZETHRESHOLD]

        # Glog transform
        features_organoids = transform_features(
            features=features_organoids, transform_type="glog", c=0.05)
        features_shrapnel = transform_features(
            features=features_shrapnel, transform_type="glog", c=0.05)

        # Load DMSO averages for the well
        # Immediately transform them as above
        dmso_avg_fn = os.path.join(
            FEATUREDIR, plate,
            "%s_dmso_summaries_%s.h5" % (plate, featuretype))
        if not os.path.isfile(dmso_avg_fn):
            calc_dmso_avg_of_plate(
                plate=plate,
                blurry_organoid_clf=blurry_organoid_clf)

        with h5py.File(dmso_avg_fn, "r") as h5handle:
            dmso_avg = transform_features(
                features=h5handle["dmso_tm05_mean_%s" % featuretype][()],
                transform_type="glog", c=0.05)
            dmso_avg_names = h5handle[
                "dmso_tm05_mean_%s_names" % featuretype][()]
            dmso_dev = transform_features(
                features=h5handle["dmso_tm05_std_%s" % featuretype][()],
                transform_type="glog", c=0.05)
            dmso_dev_names = h5handle[
                "dmso_tm05_std_%s_names" % featuretype][()]

        # Remove the number of objects entries
        dmso_avg = dmso_avg[:-1]
        dmso_dev = dmso_dev[:-1]
        dmso_avg_names = dmso_avg_names[:-1]
        dmso_dev_names = dmso_dev_names[:-1]

        features_organoids_norm = np.transpose(
            (features_organoids.transpose() - dmso_avg) /
            dmso_dev)
        features_shrapnel_norm = np.transpose(
            (features_shrapnel.transpose() - dmso_avg) /
            dmso_dev)



    return True


def cmd_learn_blurry_organoids():
    """
    This function is called when the program is called with the command line
    option "LEARN_BLURRY_ORGANOIDS"

    :return:
    """
    classifier_fn = os.path.join(
        BASEDIR, "blurry_organoid_classifier_%s.pkl" % FEATURETYPE)
    if os.path.isfile(classifier_fn):
        with open(classifier_fn, "r") as f:
            clf = pickle.load(f)
    else:
        clf = learn_blurry_organoids()
        clf["feature_type"] = FEATURETYPE
        with open(classifier_fn, "w") as f:
            pickle.dump(clf, f)
    return clf


def cmd_run_plate(plate, blurry_organoid_clf, steps_to_calc):
    """
    This function is called when the program is called with the command line
    option "RUN_PLATE"

    Go through each plate and:
        - calculate the well summaries
        - Apply glog transformation
        - Subtract the DMSO median of each plate to correct for plate effects
    The intermediate steps are stored separately to allow variations in the
    normalization to be tested and the well median calculation is a
    relatively lengthy process and unlikely to change, making it easier to
    adapt to different normalization methods (e.g. replacing glog)

    'steps_to_calc' is a tuple of steps of the calculation that should be
    performed, deleting any preexisting results out of the dataset. Recognized
    values are ("summaries", "transform", "dmso_normalization"). All other
    values are ignored. WARNING: subsequent steps expect previous steps to have
    been completed. E.g. if the expected datasets from the 'summaries' steps do
    not exist, then running only the 'transform' step will cause an error.

    :param plate:
    :param blurry_organoid_clf:
    :param steps_to_calc: A tuple
    :return:
    """
    out_fn = os.path.join(
        FEATUREDIR, plate,
        "%s_averaged_features_%s.h5" % (plate, FEATURETYPE))

    # Calc summaries
    if "summaries" in steps_to_calc:
        print("Calculating the plate summaries for '%s'" % plate)
        summary_dat = calc_plate_median(
            plate=plate, blurry_organoid_clf=blurry_organoid_clf)
        with h5py.File(out_fn, "w") as h5handle:
            h5handle.create_dataset(
                name="features_median",
                data=summary_dat["features_median"])
            h5handle.create_dataset(
                name="features_trimmed_mean_05",
                data=summary_dat["features_trimmed_mean_05"])
            h5handle.create_dataset(
                name="features_trimmed_mean_10",
                data=summary_dat["features_trimmed_mean_10"])
            h5handle.create_dataset(
                name="features_trimmed_mean_15",
                data=summary_dat["features_trimmed_mean_15"])
            h5handle.create_dataset(
                name="feature_names",
                data=summary_dat["feature_names"])
            h5handle.create_dataset(
                name="well_names",
                data=summary_dat["well_names"])

    # Add the total biomass as a feature if it doesn't exist yet
    with h5py.File(out_fn, "r") as h5handle:
        feature_names = h5handle["feature_names"][()]
    if "Total.Biomass" not in feature_names:
        with h5py.File(out_fn, "r") as h5handle:
            features_median = h5handle["features_median"][()]
            features_trimmed_05 = h5handle["features_trimmed_mean_05"][()]
            features_trimmed_10 = h5handle["features_trimmed_mean_10"][()]
            features_trimmed_15 = h5handle["features_trimmed_mean_15"][()]
            well_names = h5handle["well_names"][()]
        print("Calculating total biomass for plate '%s'" % plate)
        plate_biomass = []
        for well_name in well_names:
            seg_fn = os.path.join(
                SEGMENTATIONDIR, plate,
                "%s_DNNsegmentation.h5" % well_name)
            with h5py.File(seg_fn, "r") as h5handle:
                mask = h5handle["mask"][()]
                plate_biomass.append(np.sum(mask > 0))
        features_median = np.concatenate(
            (features_median, np.expand_dims(plate_biomass, 0)))
        features_trimmed_05 = np.concatenate(
            (features_trimmed_05, np.expand_dims(plate_biomass, 0)))
        features_trimmed_10 = np.concatenate(
            (features_trimmed_10, np.expand_dims(plate_biomass, 0)))
        features_trimmed_15 = np.concatenate(
            (features_trimmed_15, np.expand_dims(plate_biomass, 0)))
        feature_names = np.concatenate(
            (feature_names, np.array(["Total.Biomass"])))
        with h5py.File(out_fn, "r+") as h5handle:
            del h5handle["features_median"]
            del h5handle["features_trimmed_mean_05"]
            del h5handle["features_trimmed_mean_10"]
            del h5handle["features_trimmed_mean_15"]
            del h5handle["feature_names"]
            h5handle.create_dataset(
                name="features_median", data=features_median)
            h5handle.create_dataset(
                name="features_trimmed_mean_05", data=features_trimmed_05)
            h5handle.create_dataset(
                name="features_trimmed_mean_10", data=features_trimmed_10)
            h5handle.create_dataset(
                name="features_trimmed_mean_15", data=features_trimmed_15)
            h5handle.create_dataset(
                name="feature_names", data=feature_names)

    # Transform
    if "transform" in steps_to_calc:
        print("Calculating glog transform for plate '%s' ..." % plate)
        with h5py.File(out_fn, "r+") as h5handle:
            if "features_median_glog" in h5handle.keys():
                del h5handle["features_median_glog"]
            if "features_trimmed_mean_05_glog" in h5handle.keys():
                del h5handle["features_trimmed_mean_05_glog"]
            if "features_trimmed_mean_10_glog" in h5handle.keys():
                del h5handle["features_trimmed_mean_10_glog"]
            if "features_trimmed_mean_15_glog" in h5handle.keys():
                del h5handle["features_trimmed_mean_15_glog"]

            features_median = h5handle["features_median"][()]
            features_trimmed_05 = h5handle["features_trimmed_mean_05"][()]
            features_trimmed_10 = h5handle["features_trimmed_mean_10"][()]
            features_trimmed_15 = h5handle["features_trimmed_mean_15"][()]

            features_median = transform_features(
                features_median, "glog", c=0.05)
            features_trimmed_05 = transform_features(
                features_trimmed_05, "glog", c=0.05)
            features_trimmed_10 = transform_features(
                features_trimmed_10, "glog", c=0.05)
            features_trimmed_15 = transform_features(
                features_trimmed_15, "glog", c=0.05)

            h5handle.create_dataset(
                name="features_median_glog",
                data=features_median)
            h5handle.create_dataset(
                name="features_trimmed_mean_05_glog",
                data=features_trimmed_05)
            h5handle.create_dataset(
                name="features_trimmed_mean_10_glog",
                data=features_trimmed_10)
            h5handle.create_dataset(
                name="features_trimmed_mean_15_glog",
                data=features_trimmed_15)

    # Perform DMSO normalization
    if "dmso_normalization" in steps_to_calc:
        print("Calculate DMSO median normalization for plate '%s' ..." % plate)
        with h5py.File(out_fn, "r") as h5handle:
            features_median = h5handle["features_median_glog"][()]
            features_trimmed_05 = h5handle["features_trimmed_mean_05_glog"][()]
            features_trimmed_10 = h5handle["features_trimmed_mean_10_glog"][()]
            features_trimmed_15 = h5handle["features_trimmed_mean_15_glog"][()]
            well_names = h5handle["well_names"][()]
        well_ids = np.array([
            "".join(well_name.split("_")[1:3])
            for well_name in well_names])
        layout_id = plate[11:14]
        layout = pd.read_excel(
            io=os.path.join(LAYOUTDIR, "%s.xlsx" % layout_id))
        dmso_wells = layout.loc[
            layout["Product.Name"] == "DMSO",
            "Well_ID_384"].values

        m_dmso_features_median = np.nanmean(
            features_median[:, np.in1d(well_ids, dmso_wells)],
            axis=1)
        m_dmso_features_tm05 = np.nanmean(
            features_trimmed_05[:, np.in1d(well_ids, dmso_wells)],
            axis=1)
        m_dmso_features_tm10 = np.nanmean(
            features_trimmed_10[:, np.in1d(well_ids, dmso_wells)],
            axis=1)
        m_dmso_features_tm15 = np.nanmean(
            features_trimmed_15[:, np.in1d(well_ids, dmso_wells)],
            axis=1)
        s_dmso_features_median = np.nanstd(
            features_median[:, np.in1d(well_ids, dmso_wells)],
            axis=1)
        s_dmso_features_tm05 = np.nanstd(
            features_trimmed_05[:, np.in1d(well_ids, dmso_wells)],
            axis=1)
        s_dmso_features_tm10 = np.nanstd(
            features_trimmed_10[:, np.in1d(well_ids, dmso_wells)],
            axis=1)
        s_dmso_features_tm15 = np.nanstd(
            features_trimmed_15[:, np.in1d(well_ids, dmso_wells)],
            axis=1)

        features_median_norm = np.transpose(
            (features_median.transpose() - m_dmso_features_median) /
            s_dmso_features_median)
        features_tm05_norm = np.transpose(
            (features_trimmed_05.transpose() - m_dmso_features_tm05) /
            s_dmso_features_tm05)
        features_tm10_norm = np.transpose(
            (features_trimmed_10.transpose() - m_dmso_features_tm10) /
            s_dmso_features_tm10)
        features_tm15_norm = np.transpose(
            (features_trimmed_15.transpose() - m_dmso_features_tm15) /
            s_dmso_features_tm15)

        with h5py.File(out_fn, "r+") as h5handle:
            if "features_median_normalized" in h5handle.keys():
                del h5handle["features_median_normalized"]
            if "features_trimmed_mean_05_normalized" in h5handle.keys():
                del h5handle["features_trimmed_mean_05_normalized"]
            if "features_trimmed_mean_10_normalized" in h5handle.keys():
                del h5handle["features_trimmed_mean_10_normalized"]
            if "features_trimmed_mean_15_normalized" in h5handle.keys():
                del h5handle["features_trimmed_mean_15_normalized"]
            if "readme" in h5handle.keys():
                del h5handle["readme"]

            h5handle.create_dataset(
                name="features_median_normalized",
                data=features_median_norm)
            h5handle.create_dataset(
                name="features_trimmed_mean_05_normalized",
                data=features_tm05_norm)
            h5handle.create_dataset(
                name="features_trimmed_mean_10_normalized",
                data=features_tm10_norm)
            h5handle.create_dataset(
                name="features_trimmed_mean_15_normalized",
                data=features_tm15_norm)
            h5handle.create_dataset(
                name="readme", data=np.array(
                    ["Well Median -> glog -> Subtract DMSO median"]))


def cmd_run_cell_lines():
    # Go through each cell line and:
    # - Calculate the z score
    all_cell_lines = sorted(set([
        plate[0:7] for plate in os.listdir(FEATUREDIR) if
        plate.startswith("D0") or plate.startswith("M001")]))

    all_plates = sorted([
        plate for plate in os.listdir(FEATUREDIR) if
        plate.startswith("D0") or plate.startswith("M001")])

    for cell_line in all_cell_lines:
        cl_fn = os.path.join(
            BASEDIR, "%s_averaged_features_%s.h5"
                     % (cell_line, FEATURETYPE))
        if os.path.isfile(cl_fn):
            print("Cell line features file already exists for '%s'" % cell_line)
            continue

        print("Combining cell line '%s'" % cell_line)

        cl_plates = sorted([
            plate for plate in all_plates
            if plate.startswith(cell_line)])

        # Use only re-imaged plates (9XX vs 0XX)
        plate_ids = [s[8:11] for s in cl_plates]
        use_plate = []
        for plate_id in plate_ids:
            if plate_id[0] == "9":
                use_plate.append(True)
                continue
            reimaged = "9" + plate_id[1:3]
            if reimaged in plate_ids:
                use_plate.append(False)
            else:
                use_plate.append(True)
        cl_plates = [
            cl_plates[i] for i in
            range(len(cl_plates)) if use_plate[i]]

        # Sort the plates
        cl_plates = sorted(cl_plates, key=lambda x: x[9:14])

        # Set replicates
        plate_ids = [s[12:14] for s in cl_plates]
        replicates = {k: 0 for k in cl_plates}
        for plate_id in set(plate_ids):
            rep_plates = [
                cl_plates[i] for i in range(len(cl_plates))
                if plate_ids[i] == plate_id]
            rep_plates_num = [s[9:11] for s in rep_plates]
            replicates[rep_plates[rep_plates_num.index(max(rep_plates_num))]] = 2
            replicates[rep_plates[rep_plates_num.index(min(rep_plates_num))]] = 1

        # Load features and layouts
        cl_features_median = []
        cl_features_tm05 = []
        cl_features_tm10 = []
        cl_features_tm15 = []
        cl_feature_names = []
        cl_well_names = []
        cl_replicates = []
        cl_drugs = []
        cl_concentrations = []
        for plate in cl_plates:
            feature_fn = os.path.join(
                FEATUREDIR, plate,
                "%s_averaged_features_%s.h5" % (plate, FEATURETYPE))
            with h5py.File(feature_fn, "r") as h5handle:
                features_median = h5handle["features_median_normalized"][()]
                features_tm05 = h5handle["features_trimmed_mean_05_normalized"][()]
                features_tm10 = h5handle["features_trimmed_mean_10_normalized"][()]
                features_tm15 = h5handle["features_trimmed_mean_15_normalized"][()]
                feature_names = h5handle["feature_names"][()]
                well_names = h5handle["well_names"][()]
            cl_replicates.append([replicates[plate]] * len(well_names))
            cl_features_median.append(features_median)
            cl_features_tm05.append(features_tm05)
            cl_features_tm10.append(features_tm10)
            cl_features_tm15.append(features_tm15)
            cl_feature_names.append(feature_names)
            cl_well_names.append(well_names)

            layout_id = plate[11:14]
            layout = pd.read_excel(
                io=os.path.join(LAYOUTDIR, "%s.xlsx" % layout_id))
            well_ids = np.array([
                "".join(well_name.split("_")[1:3])
                for well_name in well_names])
            for well_id in well_ids:
                cl_drugs.append(layout.loc[
                                    layout["Well_ID_384"] == well_id,
                                    "Product.Name"].values)
                if "concentration" in layout.columns:
                    cl_concentrations.append(layout.loc[
                                                 layout["Well_ID_384"] == well_id,
                                                 "concentration"].values)
                else:
                    cl_concentrations.append(np.array([np.nan]))
        fname_iter = iter(cl_feature_names)
        if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
            raise Warning("Not all feature names were identical between files")
        cl_feature_names = cl_feature_names[0]
        cl_features_median = np.concatenate(cl_features_median, axis=1)
        cl_features_tm05 = np.concatenate(cl_features_tm05, axis=1)
        cl_features_tm10 = np.concatenate(cl_features_tm10, axis=1)
        cl_features_tm15 = np.concatenate(cl_features_tm15, axis=1)
        cl_well_names = np.concatenate(cl_well_names, axis=0)
        cl_replicates = np.concatenate(cl_replicates, axis=0)
        cl_drugs = np.concatenate(cl_drugs, axis=0).astype(np.str)
        cl_concentrations = np.concatenate(cl_concentrations, axis=0)

        # Calculate z score across cell line and per replicate
        feature_sets = {
            "cl_features_median": cl_features_median,
            "cl_features_tm05": cl_features_tm05,
            "cl_features_tm10": cl_features_tm10,
            "cl_features_tm15": cl_features_tm15}
        feature_sets_z = {}
        feature_sets_z_rep = {}
        for feature_set_key in feature_sets.keys():
            feature_set = feature_sets[feature_set_key]
            feature_set_masked = np.ma.array(
                feature_set, mask=np.isnan(feature_set))
            # Z score across entire cell line
            feature_sets_z[feature_set_key] = scipy.stats.mstats.zscore(
                a=feature_set_masked, axis=1)
            # Z score across replicate
            feature_set_z_rep = np.zeros(
                shape=feature_set.shape)
            feature_set_z_rep[..., cl_replicates == 1] = scipy.stats.mstats.zscore(
                    a=feature_set_masked[..., cl_replicates == 1], axis=1)
            feature_set_z_rep[..., cl_replicates == 2] = scipy.stats.mstats.zscore(
                a=feature_set_masked[..., cl_replicates == 2], axis=1)
            feature_sets_z_rep[feature_set_key] = feature_set_z_rep

        # Save data
        with h5py.File(cl_fn, "w-") as h5handle:
            h5handle.create_dataset(
                name="features_median_zscore",
                data=feature_sets_z["cl_features_median"])
            h5handle.create_dataset(
                name="features_tm05_zscore",
                data=feature_sets_z["cl_features_tm05"])
            h5handle.create_dataset(
                name="features_tm10_zscore",
                data=feature_sets_z["cl_features_tm10"])
            h5handle.create_dataset(
                name="features_tm15_zscore",
                data=feature_sets_z["cl_features_tm15"])
            h5handle.create_dataset(
                name="features_median_zscore_rep",
                data=feature_sets_z_rep["cl_features_median"])
            h5handle.create_dataset(
                name="features_tm05_zscore_rep",
                data=feature_sets_z_rep["cl_features_tm05"])
            h5handle.create_dataset(
                name="features_tm10_zscore_rep",
                data=feature_sets_z_rep["cl_features_tm10"])
            h5handle.create_dataset(
                name="features_tm15_zscore_rep",
                data=feature_sets_z_rep["cl_features_tm15"])
            h5handle.create_dataset(
                name="features_median", data=cl_features_median)
            h5handle.create_dataset(
                name="features_tm05", data=cl_features_tm05)
            h5handle.create_dataset(
                name="features_tm10", data=cl_features_tm10)
            h5handle.create_dataset(
                name="features_tm15", data=cl_features_tm15)
            h5handle.create_dataset(
                name="feature_names", data=cl_feature_names)
            h5handle.create_dataset(name="well_names", data=cl_well_names)
            h5handle.create_dataset(name="drugs", data=cl_drugs)
            h5handle.create_dataset(name="replicates", data=cl_replicates)
            h5handle.create_dataset(
                name="concentrations", data=cl_concentrations)


if __name__ == "__main__":
    # Set up command line arguments
    cmd_args = sys.argv
    if len(cmd_args) < 2:
        print(
            "Usage: %s [LEARN_BLURRY_ORGANOIDS|RUN_PLATE|RUN_CELL_LINES] "
            "<PLATE>" % cmd_args[0])
        sys.exit()

    if cmd_args[1] == "LEARN_BLURRY_ORGANOIDS":
        cmd_learn_blurry_organoids()
        sys.exit()
    elif cmd_args[1] == "RUN_PLATE":
        bo_clf = cmd_learn_blurry_organoids()
        cmd_run_plate(cmd_args[2], bo_clf, cmd_args[2:])
        sys.exit()
    elif cmd_args[1] == "RUN_CELL_LINES":
        cmd_run_cell_lines()
        pass
    else:
        print(
            "Usage: %s [LEARN_BLURRY_ORGANOIDS|RUN_PLATE"
            "|RUN_CELL_LINES] <PLATE>" % cmd_args[0])
        sys.exit()
