# This script serves as demonstration. The proper (and updated!)
# code is embedded in the function(s) for preprocessing features

import numpy as np
import sklearn.ensemble
import sklearn.model_selection
import h5py
import os
import matplotlib.pyplot as plt
import random
import re


def learn_blurry_organoids(features, feature_names):
    """
    Trains a classifier to tell the difference between blurry and in-focus
    organoids based on the features.

    The steps are:
    1. Initial annotation with a high-precision/low-recall segmentation based
       on maximal DAPI/Cy3 intensity.
    2. Train random forest classifier on all unrelated features to prevent
       overfitting

    :param features: A 2D numpy array with the shape (samples, features)
    :param feature_names: A 1D numpy array with the feature names. Must
           have the same length as the second dimension of 'features'
    :return:
        - A random forest classifier object
        - The expected accuracy on a validation set
    """

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
    # remove_features_names = [
    #     s for s in feature_names if
    #     re.match("x\.(a|ab|ac|Ba|Bab|Bac|ba|ca|aB|abB|baB|acB|caB)\.b", s)
    #     is not None]
    remove_features_names = [
        s for s in feature_names if
        re.match("x\..*\.b", s)
        is not None]
    # remove_features_names = (
    #     "x.a.b.q099", "x.a.b.q095", "x.ac.b.q099", "x.ac.b.q095",
    #     "x.ab.b.q099", "x.ab.b.q095", "x.Ba.b.q099", "x.Ba.b.q095",
    #     "x.Bac.b.q099", "x.Bac.b.q095", "x.Bab.b.q099", "x.Bab.b.q095")
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

    return rfclf, val_acc


def test_random_forest_classifier(plate, row, col, classifier, field_id=0):
    feature_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
    segmentation_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/segmentation"
    projection_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5projection"

    well_to_load = (plate, row, col)
    feature_fn = os.path.join(
        feature_dir, plate, "wells",
        "%s_%s_%s_features.h5" % well_to_load)
    with h5py.File(feature_fn, "r") as h5handle:
        features = h5handle["features_clumps"][()]
        feature_names = h5handle["feature_names_clumps"][()]

    segmask_fn = os.path.join(
        segmentation_dir, plate,
        "%s_%s_%s_DNNsegmentation.h5" % well_to_load)
    with h5py.File(segmask_fn, "r") as h5handle:
        m = h5handle["mask"][field_id, ...]
        m[m > 0] = 1
        mask = m

    proj_fn = os.path.join(
        projection_dir, plate,
        "%s_%s_%s_contrastProjections.h5" % well_to_load)
    with h5py.File(proj_fn, "r") as h5handle:
        proj = h5handle["images"][field_id, ...]

    features = features.transpose()
    features = features[features[:, -1] == (field_id + 1), ...]

    # Remove small spots as these aren't organoids and irrelevant for this
    # visualization
    # sizes = np.squeeze(test_features[:, np.where(test_feature_names == "x.0.s.area")[0]])
    # test_features = test_features[[s >= 1000 for s in sizes], ...]

    proj = np.transpose(proj, (1, 2, 0))
    proj = proj.astype(np.float32)
    proj /= np.min(np.max(proj, axis=(0, 1)))

    # remove_features_names = [
    #     s for s in feature_names if
    #     re.match("x\.(a|ab|ac|Ba|Bab|Bac|ba|ca|aB|abB|baB|acB|caB)\.b", s)
    #     is not None]
    remove_features_names = [
        s for s in feature_names if
        re.match("x\..*\.b", s)
        is not None]
    # remove_features_names = (
    #     "x.a.b.q099", "x.a.b.q095", "x.ac.b.q099", "x.ac.b.q095",
    #     "x.ab.b.q099", "x.ab.b.q095", "x.Ba.b.q099", "x.Ba.b.q095",
    #     "x.Bac.b.q099", "x.Bac.b.q095", "x.Bab.b.q099", "x.Bab.b.q095")
    index_mask = np.ma.array(range(len(feature_names)), mask=False)
    for feat_name in remove_features_names:
        feat_index = np.where(feature_names == feat_name)[0][0]
        index_mask[feat_index] = np.ma.masked

    features_masked = features[:, index_mask.compressed()]
    predicted_class = classifier.predict(features_masked)

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    plt.imshow(mask)
    for organoid_id in range(features_masked.shape[0]):
        # R and Python store x- and y- coordinates in reverse order for matrices.
        # However, matplotlib flips them AGAIN, effectively cancelling the
        # difference. The double reverse order is kept here for clarity's sake,
        # however.
        y, x = features_masked[organoid_id][0:2]
        bx, by = x, y
        # bx = np.sign(mask.shape[0]/2 - x) * 30 + x
        # by = np.sign(mask.shape[1]/2 - y) * 30 + y
        if predicted_class[organoid_id] == 0:
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


if __name__ == "__main__":
    feature_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
    segmentation_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/segmentation"
    projection_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5projection"
    blurry_well_fn = "/Users/jansauer/Thesis/Projects/PROMISE/FilterBlurryWells/" \
                     "mouse/blurry_wells_predicted.txt"

    # Choose 25 random wells from each available plate (exclude all-blurry wells)
    with open(blurry_well_fn, "r") as f:
        blurry_wells = [s.strip() for s in f.readlines()]
    all_plates = [s for s in os.listdir(feature_dir) if s.startswith("M0")]
    features = []
    feature_names = []
    for plate in all_plates:
        wells = [s for s in os.listdir(
            os.path.join(feature_dir, plate, "wells"))]
        wells = [s for s in wells if s[0:19] not in blurry_wells]
        wells = random.sample(wells, 25)
        for well in wells:
            feature_fn = os.path.join(feature_dir, plate, "wells", well)
            with h5py.File(feature_fn, "r") as h5handle:
                features.append(h5handle["features_clumps"][()])
                feature_names.append(h5handle["feature_names_clumps"][()])
    features = np.concatenate(features, axis=1)
    features = features.transpose()
    feature_names = feature_names[0]

    # Remove small spots as these aren't organoids and irrelevant for this problem
    sizes = np.squeeze(features[:, np.where(feature_names == "x.0.s.area")[0]])
    features_trimmed = features[[s >= 1000 for s in sizes], ...]

    clf, acc = learn_blurry_organoids(features_trimmed, feature_names)

    # well = ("M001W01P012L05", "N", "06")
    # well = ("M001A03P006L05", "N", "06")
    # well = ("M001B04P008L07", "M", "12")
    well = ("D018T01P906L03", "A", "12")
    test_random_forest_classifier(
        plate=well[0], row=well[1], col=well[2],
        classifier=clf, field_id=0)
