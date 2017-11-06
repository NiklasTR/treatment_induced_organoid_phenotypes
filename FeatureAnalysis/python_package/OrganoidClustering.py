import h5py
import BlurryOrganoidClassifier
import Config
import os
import numpy as np
import re
import pickle
import sys


def demo_organoid_clustering():
    """
    This is a prototype function to demo the clustering of organoids within a well
    :return:
    """
    pass


def get_organoid_labels_for_well(well_id):
    """
    Applies labels to the organoids in a well in that it tests for organoids being
    blurry and alive
    :return:
    """
    well_id_split = well_id.split("_")
    plate, row, col = well_id_split[0:3]
    features_fn = os.path.join(
        Config.FEATUREDIR, plate, "wells",
        "%s_%s_%s_features.h5" % (plate, row, col))

    # Load features
    with h5py.File(features_fn, "r") as h5handle:
        features = h5handle["features_%s" % Config.FEATURETYPE][()]
        feature_names = h5handle["feature_names_%s" % Config.FEATURETYPE][()]

    # Remove FIELD info
    if "FIELD" in feature_names:
        features = np.delete(
            features, np.where(feature_names == "FIELD")[0], axis=0)
        feature_names = np.delete(
            feature_names, np.where(feature_names == "FIELD")[0])

    # Label blurry organoids
    clf = BlurryOrganoidClassifier.get_blurry_organoid_classifier()

    features_clf = np.array(features).transpose()
    features_clf = features_clf[
        ..., [f in clf["feature_names"] for f in feature_names]]

    is_focused = clf["clf"].predict(features_clf)
    is_blurry = 1 - is_focused

    # Label organoids and shrapnel
    f_size = features[np.where(
        feature_names == "x.0.s.area")[0][0]]
    is_organoid = (f_size >= Config.SIZETHRESHOLD).astype(int)
    features_shrapnel = features[:, f_size < Config.SIZETHRESHOLD]
    features_organoids = features[:, f_size >= Config.SIZETHRESHOLD]

    # Load DMSO averages for the well
    # Immediately transform them as above
    dmso_avg_fn = os.path.join(
        Config.FEATUREDIR, plate,
        "%s_dmso_summaries_%s.h5" % (plate, Config.FEATURETYPE))
    with h5py.File(dmso_avg_fn, "r") as h5handle:
        dmso_avg_organoids = h5handle["dmso_tm05_mean_organoids"][()]
        dmso_avg_organoids_names = h5handle[
            "dmso_tm05_mean_organoids_names"][()]
        dmso_dev_organoids = h5handle["dmso_tm05_std_organoids"][()]
        dmso_dev_organoids_names = h5handle[
            "dmso_tm05_std_organoids_names"][()]
        dmso_avg_shrapnel = h5handle["dmso_tm05_mean_shrapnel"][()]
        dmso_avg_shrapnel_names = h5handle[
            "dmso_tm05_mean_shrapnel_names"][()]
        dmso_dev_shrapnel = h5handle["dmso_tm05_std_shrapnel"][()]
        dmso_dev_shrapnel_names = h5handle[
            "dmso_tm05_std_shrapnel_names"][()]

    # Remove the number of objects entries
    dmso_avg_organoids = dmso_avg_organoids[[
        re.search("num\.of\.objects", s) is None
        for s in dmso_avg_organoids_names]]
    dmso_dev_organoids = dmso_dev_organoids[[
        re.search("num\.of\.objects", s) is None
        for s in dmso_dev_organoids_names]]
    dmso_avg_shrapnel = dmso_avg_shrapnel[[
        re.search("num\.of\.objects", s) is None
        for s in dmso_avg_shrapnel_names]]
    dmso_dev_shrapnel = dmso_dev_shrapnel[[
        re.search("num\.of\.objects", s) is None
        for s in dmso_dev_shrapnel_names]]

    # Normalize features
    dmso_expected = np.stack(
        [(dmso_avg_shrapnel, dmso_avg_organoids)[i]
         for i in is_organoid], axis=1)
    dmso_variation = np.stack(
        [(dmso_dev_shrapnel, dmso_dev_organoids)[i]
         for i in is_organoid], axis=1)
    features_normalized = (features - dmso_expected) / dmso_variation

    # Classify organoids as live or dead
    cell_line = plate[0:7]
    classifier_fn = "OrganoidClassifier_%s.pkl" % cell_line
    with open(classifier_fn, "rb") as f:
        clf = pickle.load(f)

    # Reduce features to required set
    features_clf = features_normalized[np.in1d(feature_names, clf[2]), :]

    # Remove bad organoids
    bad_organoids = np.sum(np.isfinite(features_clf), axis=0)
    features_clf = features_clf[:, bad_organoids != 0]

    # Impute missing values
    features_mask = np.ma.array(features_clf, mask=False)
    features_mask[~np.isfinite(features_clf)] = np.ma.masked
    col_medians = np.nanmedian(features_mask, axis=0)
    col_medians = np.stack([col_medians] * features_clf.shape[0])
    features_clf[~np.isfinite(features_clf)] = col_medians[~np.isfinite(features_clf)]

    # Apply classifier to wells
    features_clf = np.transpose(features_clf)
    prediction = clf[0].predict(features_clf)

    label = prediction
    label[is_blurry == 1] = "BLURRY"

    return label


if __name__ == "__main__":
    cell_line = sys.argv[1]
    plates = [
        plate for plate in os.listdir(Config.FEATUREDIR)
        if plate.startswith(cell_line)]
    for plate in plates:
        wells = [
            "_".join(well.split("_")[0:3]) for well in
            os.listdir(os.path.join(Config.FEATUREDIR, plate, "wells"))
            if well.startswith(plate)]
        for well in wells:
            labels = get_organoid_labels_for_well(well)
            out_fn = os.path.join("well_labels", "%s_labels.txt" % well)
            with open(out_fn, "w") as f:
                f.writelines([label + "\n" for label in labels])
