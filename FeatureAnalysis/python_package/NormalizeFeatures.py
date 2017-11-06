# These functions normalize the organoid features
import os
import pandas as pd
import numpy as np
import re
import BlurryOrganoidClassifier
import Config
import LoadFeatures
import h5py
import Utils
import FeatureSummaries


def get_normalized_organoid_features(features, feature_names, well_names):
    """
    Applies the DMSO normalization to the features, i.e.

    f = (f - mean_dmso) / sd_dmso
    :param features:
    :param feature_names:
    :param well_names:
    :return:
    """

    wells_present = sorted(set([
        "_".join(well.split("_")[0:3])
        for well in well_names]))
    features_normalized = []
    feature_names_normalized = []
    well_names_normalized = []
    object_types = []
    for well in wells_present:
        plate = well.split("_")[0]
        normalize_features_dir = os.path.join(
            Config.FEATUREDIR, plate, "wells_normalized")
        if not os.path.isdir(normalize_features_dir):
            os.makedirs(normalize_features_dir)
        normalized_features_fn = os.path.join(
            normalize_features_dir, "%s_features_normalized.h5" % well)
        if os.path.isfile(normalized_features_fn):
            well_features = LoadFeatures.load_organoid_features(
                wells=[well], normalized=True)
            features_normalized.append(well_features["features"])
            feature_names_normalized.append(well_features["feature_names"])
            well_names_normalized.append(well_features["well_names"])
            object_types.append(well_features["object_type"])
        else:
            well_features = features[..., well_names == well]
            dmso_average = get_dmso_average_for_plate(plate)

            # Separate into organoids and shrapnel
            f_size = well_features[np.where(
                feature_names == "x.0.s.area")[0][0]]
            features_shrapnel = well_features[:, f_size < Config.SIZETHRESHOLD]
            features_organoids = well_features[:, f_size >= Config.SIZETHRESHOLD]

            # Remove the number of objects entries
            dmso_avg_organoids = dmso_average["organoids"]["expected"][[
                re.search("num\.of\.objects", s) is None
                for s in dmso_average["organoids"]["feature_names_expected"]]]
            dmso_dev_organoids = dmso_average["organoids"]["variation"][[
                re.search("num\.of\.objects", s) is None
                for s in dmso_average["organoids"]["feature_names_variation"]]]
            dmso_avg_shrapnel = dmso_average["shrapnel"]["expected"][[
                re.search("num\.of\.objects", s) is None
                for s in dmso_average["shrapnel"]["feature_names_expected"]]]
            dmso_dev_shrapnel = dmso_average["shrapnel"]["variation"][[
                re.search("num\.of\.objects", s) is None
                for s in dmso_average["shrapnel"]["feature_names_variation"]]]

            features_organoids_norm = np.transpose(
                (features_organoids.transpose() - dmso_avg_organoids) /
                dmso_dev_organoids)
            features_shrapnel_norm = np.transpose(
                (features_shrapnel.transpose() - dmso_avg_shrapnel) /
                dmso_dev_shrapnel)
            features_norm = np.concatenate(
                (features_organoids_norm, features_shrapnel_norm),
                axis=1)
            object_type = np.repeat(
                ("Organoid", "Shrapnel"),
                (features_organoids_norm.shape[1],
                 features_shrapnel_norm.shape[1]))

            features_normalized.append(features_norm)
            feature_names_normalized.append(feature_names)
            well_names_normalized.append(np.repeat([well], features_norm.shape[1]))
            object_types.append(object_type)

            with h5py.File(normalized_features_fn, "a") as h5handle:
                h5handle.create_dataset(
                    name="features_%s" % Config.FEATURETYPE,
                    data=features_norm)
                h5handle.create_dataset(
                    name="feature_names_%s" % Config.FEATURETYPE,
                    data=feature_names)
                h5handle.create_dataset(
                    name="object_type_%s" % Config.FEATURETYPE,
                    data=object_type)

    features_normalized = np.concatenate(features_normalized, axis=1)
    fname_iter = iter(feature_names_normalized)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    feature_names_normalized = feature_names_normalized[0]
    well_names_normalized = np.concatenate(well_names_normalized)
    object_types = np.concatenate(object_types)

    return {
        "features": features_normalized,
        "feature_names": feature_names_normalized,
        "well_names": well_names_normalized,
        "object_type": object_types}


def get_dmso_average_for_plate(plate):
    """
    Load the average features for the DMSO controls of a given plate
    :param plate:
    :return:
    """
    out_fn = os.path.join(
        Config.FEATUREDIR, plate,
        "%s_dmso_summaries_%s.h5" % (plate, Config.FEATURETYPE))
    if os.path.isfile(out_fn):
        dmso_averages = {}
        with h5py.File(out_fn, "r") as h5handle:
            for object_type in h5handle.keys():
                dmso_averages[object_type] = dict()
                for data_set in h5handle[object_type].keys():
                    dmso_averages[object_type][data_set] = \
                        h5handle[object_type][data_set][()]
    else:
        dmso_averages = calc_dmso_organoid_average_of_plate(plate)
        with h5py.File(out_fn, "w") as h5handle:
            for object_type in dmso_averages.keys():
                h5group = h5handle.create_group(name=object_type)
                for data_set in dmso_averages[object_type].keys():
                    h5group.create_dataset(
                        name=data_set,
                        data=dmso_averages[object_type][data_set])

    return dmso_averages


def calc_dmso_organoid_average_of_plate(plate):
    """
    Calculates the feature summaries of objects in DMSO wells on a plate on an
    object level.

    :param plate:
    :return:
    """

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

    # Remove blurry organoids
    features = BlurryOrganoidClassifier.remove_blurry_organoids(**features)

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
    organoids_in_wells = np.array(organoids_in_wells)
    shrapnel_in_wells = np.array(shrapnel_in_wells)
    # This hack creates an array of "features" with an artificial size
    # denoting whether the objects are organoids or shrapnel. This makes
    # it compatible with the function 'calc_feature_summary'
    objects_in_wells = np.stack((
        np.concatenate((organoids_in_wells, shrapnel_in_wells)),
        np.repeat(
            (Config.SIZETHRESHOLD+1, Config.SIZETHRESHOLD-1),
            (len(organoids_in_wells), len(shrapnel_in_wells)))),
        axis=0)

    # Calculate the summaries over the features and add the summaries of
    # the number of objects as a feature
    summary = FeatureSummaries.calc_feature_summary(
        features=features["features"],
        feature_names=features["feature_names"],
        summary_func_middle=Utils.trim_func,
        summary_func_var=Utils.trim_func,
        kwargs_middle={"func": np.nanmean, "percent": 0.05},
        kwargs_var={"func": np.nanstd, "percent": 0.05})
    summary_objects = FeatureSummaries.calc_feature_summary(
        features=objects_in_wells,
        feature_names=np.array(["num.of.objects", "x.0.s.area"]),
        summary_func_middle=Utils.trim_func,
        summary_func_var=Utils.trim_func,
        kwargs_middle={"func": np.nanmean, "percent": 0.05},
        kwargs_var={"func": np.nanstd, "percent": 0.05})
    for ii in ("organoids", "shrapnel"):
        for jj in ("expected", "variation"):
            summary[ii][jj] = np.append(
                summary[ii][jj],
                summary_objects[ii][jj][0])
            summary[ii]["feature_names_%s" % jj] = np.append(
                summary[ii]["feature_names_%s" % jj],
                summary_objects[ii]["feature_names_%s" % jj][0])

    return summary
