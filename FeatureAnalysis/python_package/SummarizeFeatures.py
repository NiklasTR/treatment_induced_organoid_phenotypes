import numpy as np
import statsmodels.robust
import Config
import os
import pandas as pd
import LoadFeatures
import Utils
import h5py


def calc_feature_summary(
        features, feature_names, well_names,
        summary_func_middle=np.nanmedian,
        summary_func_var=statsmodels.robust.mad,
        kwargs_middle=None,
        kwargs_var=None):
    """
    Calculate the summaries over features. Defaults to the median and MAD.
    This function splits the individual objects into organoids and shrapnel
    and treats them individually.

    :param features: A 2D numpy matrix with the shape (features, samples)
    :param feature_names:
    :param well_names:
    :param summary_func_middle:
    :param summary_func_var:
    :param kwargs_middle:
    :param kwargs_var:

    :return:
        - A feature vector for the well
        - A list of feature names
    """
    if kwargs_middle is None:
        kwargs_middle = dict()
    if kwargs_var is None:
        kwargs_var = {"center": np.nanmedian}

    # Remove the feature "FIELD"
    if "FIELD" in feature_names:
        features = np.delete(
            features, np.where(feature_names == "FIELD")[0], axis=0)
        feature_names = np.delete(
            feature_names, np.where(feature_names == "FIELD")[0])

    # Separate into organoids and shrapnel
    f_size = features[np.where(feature_names == "x.0.s.area")[0][0]]
    features_shrapnel = features[:, f_size < Config.SIZETHRESHOLD]
    features_organoids = features[:, f_size >= Config.SIZETHRESHOLD]

    # Shrapnel summaries
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
    num_shrapnel = features_shrapnel.shape[1]

    # Organoid summaries
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
    num_organoids = features_organoids.shape[1]

    return {
        "organoids": {
            "expected": features_organoids_m,
            "variation": features_organoids_v,
            "num_objects": num_organoids,
            "feature_names_expected": np.array([
                "organoids_%s_expected" % fn for fn in feature_names]),
            "feature_names_variation": np.array([
                "organoids_%s_variation" % fn for fn in feature_names])},
        "shrapnel": {
            "expected": features_shrapnel_m,
            "variation": features_shrapnel_v,
            "num_objects": num_shrapnel,
            "feature_names_expected": np.array([
                "shrapnel_%s_expected" % fn for fn in feature_names]),
            "feature_names_variation": np.array([
                "shrapnel_%s_variation" % fn for fn in feature_names])}}


def get_well_summaries(plate):
    """
    Calculate the per well summaries for a given plate.

    It calculates multiple summaries and stores them all:
    - median / mad
    - trimmed mean / trimmed sd @ 5% trimming
    - trimmed mean / trimmed sd @ 10% trimming
    - trimmed mean / trimmed sd @ 15% trimming

    :param plate:

    :return:
    """
    pass


def get_dmso_average_of_plate(plate):
    """
    Calculates the feature summaries of objects in DMSO wells on a plate.
    Treats organoids and shrapnel separately. It calculates multiple types of
    summaries: the median and the trimmed means for 5%, 10% and 15% trimming.

    If the file already exists then it is loaded from disk. If not, it is
    generated and saved.

    :param plate:
    :return:
    """

    out_fn = os.path.join(
        Config.FEATUREDIR, plate,
        "%s_dmso_summaries_%s.h5" % (plate, Config.FEATURETYPE))
    if os.path.isfile(out_fn):
        return False

    # Load the library of the plate and extract the DMSO wells names
    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))
    dmso_wells = layout.loc[
        layout["Product.Name"] == "DMSO",
        "Well_ID_384"].values

    # Load the features of the wells
    dmso_wells = [
        plate + "_" + well[0] + "_" + well[1:3]
        for well in dmso_wells]
    features = LoadFeatures.load_organoid_features(wells=dmso_wells)

    # Determine the average object numbers
    organoids_in_wells = []
    shrapnel_in_wells = []
    for dmso_well in dmso_wells:
        well_features = features["features"][
            ..., features["well_names"] == dmso_well]
        f_size = well_features[np.where(
            features["feature_names"] == "x.0.s.area")[0][0]]
        organoids_in_wells.append(np.sum(
            np.greater_equal(f_size, Config.SIZETHRESHOLD)))
        shrapnel_in_wells.append(np.sum(
            np.less(f_size, Config.SIZETHRESHOLD)))
    organoids_in_wells = np.expand_dims(organoids_in_wells, axis=0)
    shrapnel_in_wells = np.expand_dims(shrapnel_in_wells, axis=0)

    # Calculate the summaries over the features and add the summaries of
    # the number of objects as a feature
    median = calc_feature_summary(
        features=features["features"],
        feature_names=features["feature_names"],
        summary_func_middle=np.nanmedian,
        summary_func_var=statsmodels.robust.mad,
        kwargs_middle=dict(), kwargs_var={"center": np.nanmedian})
    median["organoids"]["expected"] = np.append(
        median["organoids"]["expected"],
        np.nanmedian(organoids_in_wells))
    median["organoids"]["variation"] = np.append(
        median["organoids"]["variation"],
        statsmodels.robust.mad(organoids_in_wells, axis=1))
    trimmed_05 = calc_feature_summary(
        features=features["features"],
        feature_names=features["feature_names"],
        summary_func_middle=Utils.trim_func,
        summary_func_var=Utils.trim_func,
        kwargs_middle={"func": np.mean, "percent": 0.05},
        kwargs_var={"func": np.std, "percent": 0.05})
    Utils.trim_func(organoids_in_wells, func=np.mean, percent=0.05, axis=1)
    trimmed_10 = calc_feature_summary(
        features=features["features"],
        feature_names=features["feature_names"],
        summary_func_middle=Utils.trim_func,
        summary_func_var=Utils.trim_func,
        kwargs_middle={"func": np.mean, "percent": 0.10},
        kwargs_var={"func": np.std, "percent": 0.10})
    trimmed_15 = calc_feature_summary(
        features=features["features"],
        feature_names=features["feature_names"],
        summary_func_middle=Utils.trim_func,
        summary_func_var=Utils.trim_func,
        kwargs_middle={"func": np.mean, "percent": 0.15},
        kwargs_var={"func": np.std, "percent": 0.15})

    with h5py.File(out_fn, "w") as h5handle:
        h5handle.create_dataset(
            name="dmso_median_organoids",
            data=median["organoids"]["expected"])
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