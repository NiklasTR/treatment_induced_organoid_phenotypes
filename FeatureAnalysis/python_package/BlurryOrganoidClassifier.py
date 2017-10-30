import numpy as np
import sklearn.model_selection
import sklearn.ensemble
import re
import os
import h5py
import random
import pickle
import Config


def get_blurry_organoid_classifier(classifierfn):
    """
    Calls the training function or loads an existing classifier from disk
    :return:
    """

    if os.path.isfile(classifierfn):
        with open(classifierfn, "r") as f:
            clf = pickle.load(f)
    else:
        print(
            "Training blurry organoid classifier for '%s'"
            % Config.FEATURETYPE)
        clf = learn_blurry_organoids(
            featuretype=Config.FEATURETYPE,
            blurrywellfn=Config.BLURRYWELLFN,
            featuredir=Config.FEATUREDIR)
        clf["feature_type"] = Config.FEATURETYPE
        with open(classifierfn, "w") as f:
            pickle.dump(clf, f)
    return clf


def learn_blurry_organoids(featuretype, blurrywellfn, featuredir):
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

    if featuretype not in keymap.keys():
        raise KeyError("'feature_type' must be one of '%s'"
                       % str(keymap.keys()))

    hdf5_keys = keymap[featuretype]

    with open(blurrywellfn, "r") as f:
        blurry_wells = [s.strip() for s in f.readlines()]
    all_plates = [
        s for s in os.listdir(featuredir) if
        s.startswith("M001") or s.startswith("D0")]
    features = []
    feature_names = []
    for plate in all_plates:
        wells = [s for s in os.listdir(
            os.path.join(featuredir, plate, "wells"))]
        wells = [s for s in wells if s[0:19] not in blurry_wells]
        wells = random.sample(wells, 15)
        for well in wells:
            feature_fn = os.path.join(featuredir, plate, "wells", well)
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
