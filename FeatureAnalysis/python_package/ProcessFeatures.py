# These functions process features. This includes
#   - Identifying blurry organoids
#   - Normalizing them to DMSO values
#   - Calculating well summaries

from __future__ import division
import os
import Config
import h5py
import numpy as np
import sklearn.model_selection
import sklearn.ensemble
import re
import random
import pickle
import pandas as pd
import Utils
import scipy.stats.mstats

# LOAD ORGANOIDS #


def load_organoid_features(
        wells=None, plates=None, normalized=False, strip_field_info=True):
    """
    Loads the organoid features. Can load features from all passed plates or
    all passed wells. If both 'plates' and 'wells' are not None, then 'plates'
    is ignored.

    This function is capable of loading both the normalized and raw feature
    files.

    Note that this function sorts the well and plate names, their order is NOT
    preserved.

    For strict compatibility, wells and plates must be either 'None' or a
    list/tuple, they may not be other iterable objects (e.g. string or
    dictionary)

    :param wells:
    :param plates:
    :param normalized:
    :param strip_field_info: Boolean. DO NOT CHANGE THIS unless you know
    what you're doing. The workflow currently expects this to be True.
    :return:
    """
    # Check that plates and wells are non-string iterables
    if wells is not None:
        if not isinstance(wells, (list, tuple)):
            raise ValueError("'wells' must be either None or a list/tuple")
    if plates is not None:
        if not isinstance(plates, (list, tuple)):
            raise ValueError("'plates' must be either None or a list/tuple")

    # Check feature type
    keymap = {
        "organoids": (
            "features_organoids", "feature_names_organoids",
            "object_type_organoids"),
        "clumps": (
            "features_clumps", "feature_names_clumps", "object_type_clumps")}

    if Config.FEATURETYPE not in keymap.keys():
        raise KeyError("'Config.FEATURETYPE' must be one of '%s'"
                       % str(keymap.keys()))

    hdf5_keys = keymap[Config.FEATURETYPE]

    # Set wells to load
    if wells is None and plates is None:
        raise ValueError("At least one of 'wells' and 'plates' must be a list")
    # If 'plates' is set but 'wells' is not
    elif wells is None and plates is not None:
        if normalized:
            wells = sorted([
                s for plate in plates for s in
                os.listdir(os.path.join(Config.FEATUREDIR, plate,
                                        "wells_normalized"))])
        else:
            wells = sorted([
                s for plate in plates for s in
                os.listdir(os.path.join(Config.FEATUREDIR, plate, "wells"))])
    # If 'wells' is set but 'plates' is not
    elif wells is not None and plates is None:
        # Trim any trailing text from the wells
        wells = ["_".join(well.split("_")[0:3]) for well in wells]
        if normalized:
            wells = sorted([
                s + "_features_normalized.h5" for s in wells])
        else:
            wells = sorted([s + "_features.h5" for s in wells])

    # Set well folder
    if normalized:
        well_folder = "wells_normalized"
    else:
        well_folder = "wells"

    # Load features for all wells
    features = []
    feature_names = []
    well_names = []
    object_types = []
    for well in wells:
        feature_fn = os.path.join(
            Config.FEATUREDIR, well.split("_")[0], well_folder, well)
        try:
            with h5py.File(feature_fn, "r") as h5handle:
                well_features = h5handle[hdf5_keys[0]][()]
                well_feature_names = h5handle[hdf5_keys[1]][()]
            well_names.append(np.repeat(
                "_".join(well.split("_")[0:3]),
                well_features.shape[1]))
            # well_names.append(np.repeat(
            #     re.sub("_[0-9a-zA-Z]*\.h5", "", well),
            #     well_features.shape[1]))
            features.append(well_features)
            feature_names.append(well_feature_names)
        except KeyError:
            pass

        if normalized:
            try:
                with h5py.File(feature_fn, "r") as h5handle:
                    object_types.append(h5handle[hdf5_keys[2]][()])
            except KeyError:
                pass

    # Check feature names consistency
    fname_iter = iter(feature_names)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    feature_names = feature_names[0]

    # Combine features and well names
    features = np.concatenate(features, axis=1)
    well_names = np.concatenate(well_names)

    # Remove FIELD field, this is useless for this workflow
    if strip_field_info:
        features = np.delete(
            features, np.where(feature_names == "FIELD")[0], axis=0)
        feature_names = np.delete(
            feature_names, np.where(feature_names == "FIELD")[0])

    if normalized:
        return {
            "features": features, "feature_names": feature_names,
            "well_names": well_names,
            "object_type": np.concatenate(object_types)}
    else:
        return {
            "features": features, "feature_names": feature_names,
            "well_names": well_names}


# BLURRY ORGANOIDS #


def label_blurry_organoids(features, feature_names, well_names):
    """
    Removes blurry organoids

    :param features: A 2D numpy array with the shape (features, samples)
    :param feature_names:
    :param well_names:
    :return:
    """

    clf = get_blurry_organoid_classifier()

    features_clf = np.array(features).transpose()
    features_clf = features_clf[
        ..., [f in clf["feature_names"] for f in feature_names]]

    is_focused = clf["clf"].predict(features_clf)
    object_type = np.array([
        "GOOD" if val == 1 else "BLURRY" for val in is_focused])
    # features = features[..., [s == 1 for s in is_focused]]
    # well_names = well_names[[s == 1 for s in is_focused]]

    return {
        "features": features, "feature_names": feature_names,
        "well_names": well_names, "object_type": object_type}


def get_blurry_organoid_classifier():
    """
    Loads the blurry organoid classifier
    :return:
    """

    if os.path.isfile(Config.BLURRYORGANOIDCLF):
        with open(Config.BLURRYORGANOIDCLF, "r") as f:
            clf = pickle.load(f)
    else:
        print(
            "Training blurry organoid classifier for '%s'"
            % Config.FEATURETYPE)
        clf = learn_blurry_organoids()
        clf["feature_type"] = Config.FEATURETYPE
        with open(Config.BLURRYORGANOIDCLF, "w") as f:
            pickle.dump(clf, f)
    return clf


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

    with open(Config.BLURRYWELLFN, "r") as f:
        blurry_wells = [s.strip() for s in f.readlines()]
    all_plates = [
        s for s in os.listdir(Config.FEATUREDIR) if
        s.startswith("M001") or s.startswith("D0")]
    features = []
    feature_names = []
    for plate in all_plates:
        wells = [s for s in os.listdir(
            os.path.join(Config.FEATUREDIR, plate, "wells"))]
        wells = [s for s in wells if s[0:19] not in blurry_wells]
        wells = random.sample(wells, 15)
        well_features = load_organoid_features(wells=wells)
        features.append(well_features["features"])
        feature_names.append(well_features["feature_names"])
    features = np.concatenate(features, axis=1)
    features = features.transpose()
    feature_names = feature_names[0]

    cy3_ind = np.where(feature_names == "x.a.b.q099")[0][0]
    max_cy3_intensity = features[:, cy3_ind]
    max_cy3_thresh = np.percentile(max_cy3_intensity, 95)
    min_cy3_thresh = np.percentile(max_cy3_intensity, 35)
    fitc_ind = np.where(feature_names == "x.b.b.q099")[0][0]
    max_fitc_intensity = features[:, fitc_ind]
    max_fitc_thresh = np.percentile(max_fitc_intensity, 95)
    min_fitc_thresh = np.percentile(max_fitc_intensity, 35)
    dapi_ind = np.where(feature_names == "x.c.b.q099")[0][0]
    max_dapi_intensity = features[:, dapi_ind]
    max_dapi_thresh = np.percentile(max_dapi_intensity, 95)
    min_dapi_thresh = np.percentile(max_dapi_intensity, 35)

    # The positive samples are those with pixels with intensities in the top
    # 95th percentile in at least one of the channels. Negative samples are
    # those with intensities in the bottom 35th percentile in ALL channels.
    pos_training = features[
        (max_cy3_intensity >= max_cy3_thresh) +
        (max_fitc_intensity >= max_fitc_thresh) +
        (max_dapi_intensity >= max_dapi_thresh), :]
    neg_training = features[
        (max_cy3_intensity <= min_cy3_thresh) *
        (max_dapi_intensity <= min_dapi_thresh) *
        (max_fitc_intensity <= min_fitc_thresh), :]
    min_size = min(pos_training.shape[0], neg_training.shape[0])
    pos_rand_ind = np.random.choice(
        a=range(pos_training.shape[0]), size=min_size)
    pos_training = pos_training[pos_rand_ind, :]
    neg_rand_ind = np.random.choice(
        a=range(neg_training.shape[0]), size=min_size)
    neg_training = neg_training[neg_rand_ind, :]

    x = np.concatenate((neg_training, pos_training), axis=0)
    y = np.repeat((0, 1), min_size)

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

    x = x[:, index_mask.compressed()]
    x_train, x_val, y_train, y_val = sklearn.model_selection.train_test_split(
        x, y, test_size=0.25)

    rfclf = sklearn.ensemble.RandomForestClassifier(n_estimators=100)
    rfclf.fit(x_train, y_train)
    val_acc = rfclf.score(x_val, y_val)

    return {
        "clf": rfclf, "feature_names": altered_feature_names,
        "accuracy": val_acc}


def test_blurry_organoid_classifier(well_id=None, field_id=0):
    """
    A function to visually test the blurry organoid classifier
    
    Wells to visualize:
    M001W01P012L05_N_06
    M001A03P006L05_N_06
    M001B04P008L07_M_12
    D018T01P906L03_A_12

    :param well_id:
    :param field_id:
    :return:
    """
    import matplotlib.pyplot as plt

    if well_id is None:
        well_id = "M001W01P012L05_N_06"

    plate, row, col = well_id.split("_")[0:3]

    segmentation_dir = os.path.join(
        os.path.dirname(Config.FEATUREDIR), "segmentation")
    projection_dir = os.path.join(
        os.path.dirname(Config.FEATUREDIR), "hdf5projection")

    features = load_organoid_features(wells=[well_id], strip_field_info=False)

    segmask_fn = os.path.join(
        segmentation_dir, plate,
        "%s_%s_%s_DNNsegmentation.h5" % (plate, row, col))
    with h5py.File(segmask_fn, "r") as h5handle:
        m = h5handle["mask"][field_id, ...]
        m[m > 0] = 1
        mask = m

    proj_fn = os.path.join(
        projection_dir, plate,
        "%s_%s_%s_contrastProjections.h5" % (plate, row, col))
    with h5py.File(proj_fn, "r") as h5handle:
        proj = h5handle["images"][field_id, ...]

    well_objects = features["features"][
        features["feature_names"] == "FIELD", :][0] == (field_id + 1)
    features["features"] = features["features"][..., well_objects]
    features["well_names"] = features["well_names"][well_objects]

    # Remove small spots as these aren't organoids and irrelevant for this
    # visualization
    # sizes = np.squeeze(test_features[:, np.where(test_feature_names == "x.0.s.area")[0]])
    # test_features = test_features[[s >= 1000 for s in sizes], ...]

    # Run classifier
    features = label_blurry_organoids(**features)

    # Preprocess projection for visualization
    proj = np.transpose(proj, (1, 2, 0))
    proj = proj.astype(np.float32)
    proj /= np.min(np.max(proj, axis=(0, 1)))

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    plt.imshow(mask)
    for organoid_id in range(features["features"].shape[1]):
        # R and Python store x- and y- coordinates in reverse order for matrices.
        # However, matplotlib flips them AGAIN, effectively cancelling the
        # difference. The double reverse order is kept here for clarity's sake,
        # however.
        y, x = features["features"][..., organoid_id][0:2]
        bx, by = x, y
        # bx = np.sign(mask.shape[0]/2 - x) * 30 + x
        # by = np.sign(mask.shape[1]/2 - y) * 30 + y
        if features["object_type"][organoid_id] == "BLURRY":
            # ax1.annotate("B" + str(organoid_id), xy=(y, x))
            bbox_props = dict(
                boxstyle="round", fc="r", ec="0.6", alpha=0.3)
            ax1.text(
                by, bx, "B", ha="center",
                va="center", size=10, bbox=bbox_props)
        else:
            # ax1.annotate("S" + str(organoid_id), xy=(y, x))
            bbox_props = dict(
                boxstyle="round", fc="g", ec="0.6", alpha=0.3)
            ax1.text(
                by, bx, "S", ha="center",
                va="center", size=10, bbox=bbox_props)
    ax2 = fig.add_subplot(122)
    plt.imshow(proj)
    plt.show()


def create_blurry_organoid_statistics(plate):
    """
    Creates a dataset indicating how many blurry objects there are in
    each well. Handles shrapnel and organoids individually.
    :param plate:
    :return:
    """

    out_dir = os.path.join(Config.BASEDIR, "blurry_organoid_statistics")
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    out_fn = os.path.join(out_dir, "%s_blurry_organoid_statistics.csv" % plate)
    if os.path.isfile(out_fn):
        return False

    clf = get_blurry_organoid_classifier()

    wells = sorted([
        well for well in os.listdir(os.path.join(
            Config.FEATUREDIR, plate, "wells_normalized"))])
    plate_stats = []
    well_ids = []
    plate_stats_colnames = (
        "Focused_Organoids", "Blurry_Organoids",
        "Focused_Shrapnel", "Blurry_Shrapnel")
    for well in wells:
        well_id = "_".join(well.split("_")[0:3])
        features = load_organoid_features(wells=[well_id])
        features = label_blurry_organoids(**features)
        f_size = features["features"][np.where(
            features["feature_names"] == "x.0.s.area")[0][0]]
        organoid_type = np.array(["Shrapnel"] * len(f_size))
        organoid_type[f_size >= Config.SIZETHRESHOLD] = "Organoid"
        focused_organoids = np.sum(np.logical_and(
            features["object_type"] == "GOOD",
            organoid_type == "Organoid"))
        blurry_organoids = np.sum(np.logical_and(
            features["object_type"] == "BLURRY",
            organoid_type == "Organoid"))
        focused_shrapnel = np.sum(np.logical_and(
            features["object_type"] == "GOOD",
            organoid_type == "Shrapnel"))
        blurry_shrapnel = np.sum(np.logical_and(
            features["object_type"] == "BLURRY",
            organoid_type == "Shrapnel"))
        plate_stats.append((
            focused_organoids, blurry_organoids,
            focused_shrapnel, blurry_shrapnel))
        well_ids.append(well_id)

    plate_stats = pd.DataFrame(
        plate_stats, index=well_ids, columns=plate_stats_colnames)
    plate_stats.to_csv(out_fn)

    return True


# NORMALIZE FEATURES TO DMSO AVERAGE #


def calc_feature_summary(
        features, feature_names,
        summary_func_middle,
        summary_func_var,
        kwargs_middle,
        kwargs_var,
        size_threshold=Config.SIZETHRESHOLD):
    """
    Calculate the summaries over features. Defaults to the median and MAD.
    This function splits the individual objects into organoids and shrapnel
    and treats them individually.

    :param features: A 2D numpy matrix with the shape (features, samples)
    :param feature_names:
    :param summary_func_middle:
    :param summary_func_var:
    :param kwargs_middle:
    :param kwargs_var:
    :param size_threshold:

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
    features_shrapnel = features[:, f_size < size_threshold]
    features_organoids = features[:, f_size >= size_threshold]

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


def get_normalized_organoid_features(features, feature_names, well_names, object_type):
    """
    Applies the DMSO normalization to the features, i.e.

    f = (f - mean_dmso) / sd_dmso
    :param features:
    :param feature_names:
    :param well_names:
    :param object_type:
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
            well_features = load_organoid_features(
                wells=[well], normalized=True)
            features_normalized.append(well_features["features"])
            feature_names_normalized.append(well_features["feature_names"])
            well_names_normalized.append(well_features["well_names"])
            object_types.append(well_features["object_type"])
        else:
            well_features = features[..., well_names == well]
            well_object_types = object_type[well_names == well]

            # Load DMSO averages
            dmso_average = get_dmso_average_for_plate(plate)

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

            # Build a normalization matrix fitting to the organoid type and
            # normalize the features
            f_size = well_features[np.where(
                feature_names == "x.0.s.area")[0][0]]
            object_type = np.stack([
                "Organoid" if fs >= Config.SIZETHRESHOLD else
                "Shrapnel" for fs in f_size])
            object_type[well_object_types == "BLURRY"] = "BLURRY"
            dmso_avg_matrix = np.stack([
                dmso_avg_organoids if fs >= Config.SIZETHRESHOLD else
                dmso_avg_shrapnel for fs in f_size])
            dmso_dev_matrix = np.stack([
                dmso_dev_organoids if fs >= Config.SIZETHRESHOLD else
                dmso_dev_shrapnel for fs in f_size])
            features_norm = np.transpose(
                (well_features.transpose() - dmso_avg_matrix) /
                dmso_dev_matrix)

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
    features = load_organoid_features(wells=dmso_wells)

    # Remove blurry organoids
    features = label_blurry_organoids(**features)
    focused_organoids = np.array([
        val == "GOOD" for val in features["object_type"]])
    features["features"] = features["features"][
        ..., focused_organoids]
    features["well_names"] = features["well_names"][focused_organoids]

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
    summary = calc_feature_summary(
        features=features["features"],
        feature_names=features["feature_names"],
        summary_func_middle=Utils.trim_func,
        summary_func_var=Utils.trim_func,
        kwargs_middle={"func": np.nanmean, "percent": 0.05},
        kwargs_var={"func": np.nanstd, "percent": 0.05})
    summary_objects = calc_feature_summary(
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


# SUMMARIZE WELL FEATURES #


def calc_well_summaries(plate):
    """
    Average the normalized organoid features on a plate separately for
    organoids and shrapnel

    :param plate:

    :return:
    """
    out_fn = os.path.join(
        Config.FEATUREDIR, plate,
        "%s_averaged_features_%s.h5" % (plate, Config.FEATURETYPE))
    if os.path.isfile(out_fn):
        return False

    wells = sorted([
        well for well in os.listdir(os.path.join(
            Config.FEATUREDIR, plate, "wells_normalized"))])

    features_combined_all = []
    features_combined_all_names = []
    for well in wells:
        well_id = "_".join(well.split("_")[0:3])
        features = load_organoid_features(
            wells=[well_id], normalized=True)

        features_organoids = features["features"][
            ..., features["object_type"] == "Organoid"]
        features_shrapnel = features["features"][
            ..., features["object_type"] == "Shrapnel"]

        features_organoids_summary = calc_feature_summary(
            features=features_organoids,
            feature_names=features["feature_names"],
            summary_func_middle=Utils.trim_func,
            summary_func_var=Utils.trim_func,
            kwargs_middle={"func": np.nanmean, "percent": 0.05},
            kwargs_var={"func": np.nanstd, "percent": 0.05},
            size_threshold=-np.inf)["organoids"]

        features_shrapnel_summary = calc_feature_summary(
            features=features_shrapnel,
            feature_names=features["feature_names"],
            summary_func_middle=Utils.trim_func,
            summary_func_var=Utils.trim_func,
            kwargs_middle={"func": np.nanmean, "percent": 0.05},
            kwargs_var={"func": np.nanstd, "percent": 0.05},
            size_threshold=np.inf)["shrapnel"]

        # Load DMSO averages for the number of objects
        dmso_summary = get_dmso_average_for_plate(plate)
        dmso_mean_num_organoids = dmso_summary["organoids"]["expected"][
            dmso_summary["organoids"]["feature_names_expected"] ==
            "organoids_num.of.objects_expected"]
        dmso_std_num_organoids = dmso_summary["organoids"]["variation"][
            dmso_summary["organoids"]["feature_names_variation"] ==
            "organoids_num.of.objects_variation"]
        dmso_mean_num_shrapnel = dmso_summary["shrapnel"]["expected"][
            dmso_summary["shrapnel"]["feature_names_expected"] ==
            "shrapnel_num.of.objects_expected"]
        dmso_std_num_shrapnel = dmso_summary["shrapnel"]["variation"][
            dmso_summary["shrapnel"]["feature_names_variation"] ==
            "shrapnel_num.of.objects_variation"]

        # Combine the expected and variation into a single array
        features_organoids_combined = np.concatenate((
            features_organoids_summary["expected"],
            features_organoids_summary["variation"],
            (features_organoids_summary["num_objects"] -
             dmso_mean_num_organoids) / dmso_std_num_organoids))
        features_shrapnel_combined = np.concatenate((
            features_shrapnel_summary["expected"],
            features_shrapnel_summary["variation"],
            (features_shrapnel_summary["num_objects"] -
             dmso_mean_num_shrapnel) / dmso_std_num_shrapnel))
        features_combined = np.concatenate((
            features_organoids_combined,
            features_shrapnel_combined))

        features_organoids_combined_names = np.concatenate((
            features_organoids_summary["feature_names_expected"],
            features_organoids_summary["feature_names_variation"],
            np.array(["organoids_num.of.objects"])))
        features_shrapnel_combined_names = np.concatenate((
            features_shrapnel_summary["feature_names_expected"],
            features_shrapnel_summary["feature_names_variation"],
            np.array(["shrapnel_num.of.objects"])))
        features_combined_names = np.concatenate((
            features_organoids_combined_names,
            features_shrapnel_combined_names))

        # Calculate Biomass in well
        seg_fn = os.path.join(
            Config.SEGMENTATIONDIR, plate,
            "%s_DNNsegmentation.h5" % well_id)
        with h5py.File(seg_fn, "r") as h5handle:
            mask = h5handle["mask"][()]
            well_biomass = np.sum(mask > 0)
        features_combined = np.append(features_combined, well_biomass)
        features_combined_names = np.append(
            features_combined_names, np.array(["Total.Biomass"]))

        features_combined_all.append(features_combined)
        features_combined_all_names.append(features_combined_names)

    features_combined_all = np.stack(features_combined_all, axis=1)
    fname_iter = iter(features_combined_all_names)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    features_combined_all_names = features_combined_all_names[0]

    wells = ["_".join(well.split("_")[0:3]) for well in wells]

    # Normalize the biomass by the DMSO summaries
    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))
    dmso_wells = layout.loc[
        layout["Product.Name"] == "DMSO",
        "Well_ID_384"].values
    dmso_wells = [
        "%s_%s_%s" % (plate, str(dmso_well[0]), str(dmso_well[1:3]))
        for dmso_well in dmso_wells]
    dmso_biomass = features_combined_all[
        features_combined_all_names == "Total.Biomass",
        np.in1d(wells, dmso_wells)]
    dmso_biomass_expected = Utils.trim_func(
        a=np.expand_dims(dmso_biomass, axis=0),
        func=np.nanmean, axis=1, percent=0.05)
    dmso_biomass_variation = Utils.trim_func(
        a=np.expand_dims(dmso_biomass, axis=0),
        func=np.nanstd, axis=1, percent=0.05)
    adj_biomass = (
        features_combined_all[features_combined_all_names == "Total.Biomass"] -
        dmso_biomass_expected) / dmso_biomass_variation
    features_combined_all[
        features_combined_all_names == "Total.Biomass",
        ...] = adj_biomass

    with h5py.File(out_fn, "w") as h5handle:
        h5handle.create_dataset(
            name="features_%s" % Config.FEATURETYPE,
            data=features_combined_all)
        h5handle.create_dataset(
            name="feature_names_%s" % Config.FEATURETYPE,
            data=features_combined_all_names)
        h5handle.create_dataset(
            name="well_names",
            data=wells)

    return True


def calc_cell_line_features(cell_line):
    """
    Combine and z-Normalize the features for a given cell line
    Z-scoring is done on a single replicate over a single cell line
    :param cell_line:
    :return:
    """

    cl_dir = os.path.join(Config.BASEDIR, "cell_line_features")
    if not os.path.isdir(cl_dir):
        os.makedirs(cl_dir)

    cl_fn = os.path.join(
        cl_dir, "%s_averaged_features_%s.h5"
                % (cell_line, Config.FEATURETYPE))
    if os.path.isfile(cl_fn):
        print("Cell line features file already exists for '%s'" % cell_line)
        return False

    cl_plates = sorted([
        plate for plate in os.listdir(Config.FEATUREDIR)
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
    cl_features = []
    cl_feature_names = []
    cl_well_names = []
    cl_replicates = []
    cl_drugs = []
    cl_concentrations = []
    for plate in cl_plates:
        feature_fn = os.path.join(
            Config.FEATUREDIR, plate,
            "%s_averaged_features_%s.h5" % (plate, Config.FEATURETYPE))
        with h5py.File(feature_fn, "r") as h5handle:
            features = h5handle["features_%s" % Config.FEATURETYPE][()]
            feature_names = h5handle["feature_names_%s" % Config.FEATURETYPE][()]
            well_names = h5handle["well_names"][()]
        cl_replicates.append([replicates[plate]] * len(well_names))
        cl_features.append(features)
        cl_feature_names.append(feature_names)
        cl_well_names.append(well_names)

        layout_id = plate[11:14]
        layout = pd.read_excel(
            io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))
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
    cl_features = np.concatenate(cl_features, axis=1)
    cl_well_names = np.concatenate(cl_well_names, axis=0)
    cl_replicates = np.concatenate(cl_replicates, axis=0)
    cl_drugs = np.concatenate(cl_drugs, axis=0).astype(np.str)
    cl_concentrations = np.concatenate(cl_concentrations, axis=0)

    # Calculate z score across cell line and per replicate
    features_masked = np.ma.array(
        cl_features, mask=np.isnan(cl_features))
    # Z score across entire cell line
    features_z = scipy.stats.mstats.zscore(
        a=features_masked, axis=1)
    # Z score across replicate
    features_z_rep = np.zeros(
        shape=cl_features.shape)
    features_z_rep[..., cl_replicates == 1] = scipy.stats.mstats.zscore(
        a=features_masked[..., cl_replicates == 1], axis=1)
    features_z_rep[..., cl_replicates == 2] = scipy.stats.mstats.zscore(
        a=features_masked[..., cl_replicates == 2], axis=1)
    features_z_rep = features_z_rep

    # Save data
    with h5py.File(cl_fn, "w-") as h5handle:
        h5handle.create_dataset(
            name="features_zscore_%s" % Config.FEATURETYPE,
            data=features_z)
        h5handle.create_dataset(
            name="features_zscore_rep_%s" % Config.FEATURETYPE,
            data=features_z_rep)
        h5handle.create_dataset(
            name="features_%s" % Config.FEATURETYPE,
            data=cl_features)
        h5handle.create_dataset(
            name="feature_names",
            data=cl_feature_names)
        h5handle.create_dataset(name="well_names", data=cl_well_names)
        h5handle.create_dataset(name="drugs", data=cl_drugs)
        h5handle.create_dataset(name="replicates", data=cl_replicates)
        h5handle.create_dataset(
            name="concentrations", data=cl_concentrations)
