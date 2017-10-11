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


def load_cell_line_features(
        cell_line, feature_dir, feature_type="clumps"):
    """
    Load all features of a cell line into a single array
    :param cell_line:
    :param feature_dir:
    :param feature_type
    :return:
    """

    plates = sorted([
        s for s in os.listdir(feature_dir)
        if s.startswith(cell_line)])

    # Use only re-imaged plates (9XX vs 0XX)
    plate_ids = [s[8:11] for s in plates]
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
    plates = [plates[i] for i in range(len(plates)) if use_plate[i]]

    # Load features for plates
    features = []
    feature_names = []
    well_names = []
    for plate in plates:
        feature_fn = os.path.join(
            feature_dir, plate, "%s_features.h5" % plate)
        with h5py.File(feature_fn, "r") as h5handle:
            features.append(
                h5handle["features_%s" % feature_type][()])
            feature_names.append(
                h5handle["feature_names_%s" % feature_type][()])
            well_names.append(
                [plate + "_" + wn for wn in
                 h5handle["well_names_%s" % feature_type][()]])

    fname_iter = iter(feature_names)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    feature_names = feature_names[0]
    well_names = np.concatenate(well_names)
    features = np.concatenate(features, axis=1)

    return features, feature_names, well_names


def trim_func(func, a, percent, axis):
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


def learn_blurry_organoids(blurry_well_fn, feature_dir, feature_type):
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

    :param blurry_well_fn:
    :param feature_dir:
    :param feature_type:
    :return:
        - A random forest classifier object
        - The names of features used to train the RF
        - The expected accuracy on a validation set
    """

    keymap = {
        "organoids": ("features", "feature_names"),
        "clumps": ("features_clumps", "feature_names_clumps")}

    if feature_type not in keymap.keys():
        raise KeyError("'feature_type' must be one of '%s'"
                       % str(keymap.keys()))

    hdf5_keys = keymap[feature_type]

    with open(blurry_well_fn, "r") as f:
        blurry_wells = [s.strip() for s in f.readlines()]
    all_plates = [
        s for s in os.listdir(feature_dir) if
        s.startswith("M001") or s.startswith("D0")]
    features = []
    feature_names = []
    for plate in all_plates:
        wells = [s for s in os.listdir(
            os.path.join(feature_dir, plate, "wells"))]
        wells = [s for s in wells if s[0:19] not in blurry_wells]
        wells = random.sample(wells, 15)
        for well in wells:
            feature_fn = os.path.join(feature_dir, plate, "wells", well)
            with h5py.File(feature_fn, "r") as h5handle:
                features.append(h5handle[hdf5_keys[0]][()])
                feature_names.append(h5handle[hdf5_keys[1]][()])
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


def calc_well_average(features, feature_names, blurry_organoid_clf):
    """
    Calculate the trimmed mean and standard deviation of a well's organoid
    features.

    'blurry_well_clf' should be the dictionary object returned by the
    function 'learn_blurry_organoids()'

    :param features: A 2D numpy matrix with the shape (features, samples)
    :param blurry_organoid_clf: A dictionary.
    :return:
        - A feature vector for the well
        - A list of feature names
    """

    # These are the summary functions. They could later be turned into
    # function parameters
    middle_func = np.median
    var_func = statsmodels.robust.mad

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

    # Separate into organoids and shrapnel.
    size_threshold = 2500
    f_size = features[np.where(feature_names == "x.0.s.area")[0][0]]
    features_shrapnel = features[:, f_size < size_threshold]
    features_organoids = features[:, f_size >= size_threshold]

    # Calculate feature summaries
    if features_shrapnel.shape[1] == 0:
        features_shrapnel_m = np.zeros(
            shape=features_shrapnel.shape[0])
        features_shrapnel_v = np.zeros(
            shape=features_shrapnel.shape[0])
    else:
        features_shrapnel_m = middle_func(
            features_shrapnel, axis=1)
        features_shrapnel_v = var_func(
            features_shrapnel, axis=1)
    features_shrapnel_summary = np.concatenate(
        (features_shrapnel_m, features_shrapnel_v,
         (features_shrapnel.shape[1],)))
    
    if features_organoids.shape[1] == 0:
        features_organoids_m = np.zeros(
            shape=features_organoids.shape[0])
        features_organoids_v = np.zeros(
            shape=features_organoids.shape[0])
    else:
        features_organoids_m = middle_func(
            features_organoids, axis=1)
        features_organoids_v = var_func(
            features_organoids, axis=1)
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
        if "c" not in kwargs.keys():
            c = 0.05
        return np.log((features + np.sqrt(features ** 2 + c ** 2)) / 2.0)
    else:
        raise ValueError("Unsupported transform_type ('%s')" % transform_type)


if __name__ == "__main__":
    basedir = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/preprocessing"
    feature_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
    blurry_well_fn = os.path.join(basedir, "blurry_wells_predicted.txt")
    feature_type = "clumps"

    # Begin by loading or learning a classifier to filter out blurry wells
    classifier_fn = os.path.join(
        basedir, "blurry_organoid_classifier_%s.pkl" % feature_type)
    if os.path.isfile(classifier_fn):
        with open(classifier_fn, "r") as f:
            blurry_organoid_clf = pickle.load(f)
    else:
        print("Training blurry organoid classifier")
        blurry_organoid_clf = learn_blurry_organoids(
            blurry_well_fn=blurry_well_fn,
            feature_dir=feature_dir,
            feature_type=feature_type)
        blurry_organoid_clf["feature_type"] = feature_type
        with open(classifier_fn, "w") as f:
            pickle.dump(blurry_organoid_clf, f)

    # Go through each plate and calculate the well averages
    all_plates = [
        plate for plate in os.listdir(feature_dir) if
        plate.startswith("D0") or plate.startswith("M001")]

    for plate in all_plates:
        out_fn = os.path.join(
            feature_dir, plate,
            "%s_averaged_features_%s.h5" % (plate, feature_type))
        if os.path.isfile(out_fn):
            print("Plate '%s' already processed" % plate)
            continue
        print("Processing plate '%s' ..." % plate)

        wells = sorted([
            well for well in
            os.listdir(os.path.join(feature_dir, plate, "wells"))
            if well.startswith(plate) if well.endswith(".h5")])
        avg_features = []
        avg_feature_names = []
        for well in wells:
            well_fn = os.path.join(feature_dir, plate, "wells", well)
            with h5py.File(well_fn, "r") as h5handle:
                features = h5handle["features_%s" % feature_type][()]
                feature_names = h5handle["feature_names_%s" % feature_type][()]
            avg = calc_well_average(
                features=features, feature_names=feature_names,
                blurry_organoid_clf=blurry_organoid_clf)
            avg_features.append(avg[0])
            avg_feature_names.append(avg[1])

        fname_iter = iter(avg_feature_names)
        if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
            raise Warning("Not all feature names were identical between files")
        avg_feature_names = avg_feature_names[0]
        avg_features = np.stack(avg_features, axis=1)

        well_names = np.array([well[0:19] for well in wells])

        with h5py.File(out_fn, "w") as h5handle:
            h5handle.create_dataset(name="features", data=avg_features)
            h5handle.create_dataset(name="feature_names", data=avg_feature_names)
            h5handle.create_dataset(name="well_names", data=well_names)
