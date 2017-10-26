# This script classifies individual organoids from negative and positive
# controls

from __future__ import division
import os
import h5py
import numpy as np
import sklearn.model_selection
import sklearn.ensemble
import pickle
import pandas as pd
import sys

BASEDIR = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/organoid_classifier"
FEATUREDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
LAYOUTDIR = "/collab-ag-fischer/PROMISE/layouts/python_friendly"
# SEGMENTATIONDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/segmentation"
# BLURRYWELLFN = os.path.join(BASEDIR, "blurry_wells_predicted.txt")
FEATURETYPE = "organoids"
SIZETHRESHOLD = 2500


def get_wells_for_treatment(plate, treatments=None, concentrations=None):
    """
    Get the well IDs corresponding to a given treatment
    :param plate:
    :param treatments:
    :param concentrations:
    :return:
    """

    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(LAYOUTDIR, "%s.xlsx" % layout_id))
    treatment_indices = np.repeat([True], len(layout))
    if treatments is not None:
        treatment_indices = np.in1d(
           layout["Product.Name"].values, treatments)
    concentration_indices = np.repeat([True], len(layout))
    if concentrations is not None:
        # In case the layout doesn't have a concentration column
        try:
            concentration_indices = np.in1d(
                layout["concentration"].values, concentrations)
        except KeyError:
            pass
    indices = np.logical_and(concentration_indices, treatment_indices)

    return layout.loc[indices, "Well_ID_384"].values.astype(str)


def train_classifier(cell_line, return_data=False, save=True):
    """
    Train a classifier to learn the difference between a live and a dead
    organoid. The live and dead organoids are taken from "control" wells.

    Positive controls are, as determined visually and via initial feature
    analysis, Bortezomib and Irinotecan / SN-38 at concentrations of 0.2
    and 1.0. Negative controls are DMSO wells.

    :param cell_line:
    :param return_data: Boolean. Whether or not to return the training
        data as well.
    :param save: Boolean. If True, the function tries to load a preexisting
        classifier and also saves a newly trained classifier to disk.
    :return: A tuple (classifier, accuracy, features used) or
        (classifier, accuracy, features used, data_x, data_y) if
        return_data == True
    """
    if len(cell_line) != 7:
        raise ValueError("'cell_line' should be a 7-char string")

    clf_fn = os.path.join(BASEDIR, "OrganoidClassifier_%s.pkl" % cell_line)
    if save and os.path.isfile(clf_fn):
        with open(clf_fn, "rb") as f:
            return pickle.load(f)

    # Define controls
    pos_ctrls = [
        ("Bortezomib", (0.2, 1.0)),
        ("Irinotecan / SN-38", (0.2, 1.0))]
    neg_ctrls = [("DMSO", None)]

    # Run over one cell line
    plates = [
        plate for plate in os.listdir(FEATUREDIR)
        if plate.startswith(cell_line) if plate.endswith("L08")]

    # Load control features on plate
    pos_ctrl_features = []
    pos_ctrl_feature_names = []
    pos_ctrl_object_types = []
    neg_ctrl_features = []
    neg_ctrl_feature_names = []
    neg_ctrl_object_types = []
    for plate in plates:
        pos_ctrl_wells = []
        for pos_ctrl in pos_ctrls:
            pos_ctrl_wells.append(get_wells_for_treatment(
                plate=plate, treatments=pos_ctrl[0],
                concentrations=pos_ctrl[1]))
        pos_ctrl_wells = np.concatenate(pos_ctrl_wells)
        if len(pos_ctrl_wells) == 0:
            continue

        neg_ctrl_wells = []
        for neg_ctrl in neg_ctrls:
            neg_ctrl_wells.append(get_wells_for_treatment(
                plate=plate, treatments=neg_ctrl[0],
                concentrations=neg_ctrl[1]))
        neg_ctrl_wells = np.concatenate(neg_ctrl_wells)
        if len(neg_ctrl_wells) == 0:
            continue

        for well in pos_ctrl_wells:
            features_fn = os.path.join(
                FEATUREDIR, plate, "wells_normalized",
                "%s_%s_%s_features_normalized.h5" %
                (plate, well[0], well[1:3]))
            with h5py.File(features_fn, "r") as h5handle:
                pos_ctrl_features.append(
                    h5handle["features_%s" % FEATURETYPE][()])
                pos_ctrl_feature_names.append(
                    h5handle["feature_names_%s" % FEATURETYPE][()])
                pos_ctrl_object_types.append(
                    h5handle["object_type_%s" % FEATURETYPE][()])

        for well in neg_ctrl_wells:
            features_fn = os.path.join(
                FEATUREDIR, plate, "wells_normalized",
                "%s_%s_%s_features_normalized.h5" %
                (plate, well[0], well[1:3]))
            with h5py.File(features_fn, "r") as h5handle:
                neg_ctrl_features.append(
                    h5handle["features_%s" % FEATURETYPE][()])
                neg_ctrl_feature_names.append(
                    h5handle["feature_names_%s" % FEATURETYPE][()])
                neg_ctrl_object_types.append(
                    h5handle["object_type_%s" % FEATURETYPE][()])

    # Check feature names
    pos_iter = iter(pos_ctrl_feature_names + neg_ctrl_feature_names)
    if not all(np.array_equal(next(pos_iter), rest) for rest in pos_iter):
        raise Warning("Not all feature names were identical between files")
    feature_names = pos_ctrl_feature_names[0]

    pos_ctrl_features = np.transpose(
        np.concatenate(pos_ctrl_features, axis=1))
    neg_ctrl_features = np.transpose(
        np.concatenate(neg_ctrl_features, axis=1))

    # Remove bad organoids from dataset (i.e. no features at all)
    bad_wells_pos = np.sum(np.isfinite(pos_ctrl_features), axis=1)
    pos_ctrl_features = pos_ctrl_features[bad_wells_pos >= 0, :]
    bad_wells_neg = np.sum(np.isfinite(neg_ctrl_features), axis=1)
    neg_ctrl_features = neg_ctrl_features[bad_wells_neg >= 0, :]

    # Remove bad features (any NaN values)
    combined_posneg = np.concatenate(
        (pos_ctrl_features, neg_ctrl_features), axis=0)
    bad_features = np.sum(~np.isfinite(combined_posneg), axis=0)
    pos_ctrl_features = pos_ctrl_features[:, bad_features == 0]
    neg_ctrl_features = neg_ctrl_features[:, bad_features == 0]
    feature_names = feature_names[bad_features == 0]

    # Subsample to equalize training groups
    num_samples = min(pos_ctrl_features.shape[0], neg_ctrl_features.shape[0])
    pos_ctrl_features = pos_ctrl_features[np.random.choice(
        a=np.arange(pos_ctrl_features.shape[0]), size=num_samples,
        replace=False), :]
    neg_ctrl_features = neg_ctrl_features[np.random.choice(
        a=np.arange(neg_ctrl_features.shape[0]), size=num_samples,
        replace=False), :]

    # Create training data
    x = np.concatenate((pos_ctrl_features, neg_ctrl_features), axis=0)
    y = np.repeat(("POS", "NEG"), num_samples)
    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(
        x, y, test_size=0.25)

    # Run classifier
    clf = sklearn.ensemble.RandomForestClassifier(n_estimators=100)
    clf.fit(x_train, y_train)
    acc = clf.score(x_test, y_test)

    if save:
        with open(clf_fn, "wb") as f:
            pickle.dump((clf, acc, feature_names), f, pickle.HIGHEST_PROTOCOL)

    if return_data:
        return clf, acc, feature_names, x, y
    else:
        return clf, acc, feature_names


@DeprecationWarning
def select_minimal_features_for_classifier(
        clf, x, y, feature_names, nseq=None):
    """
    Retrain a classifier with only the most important features and determine
    the validation accuracies for each subset of features
    :param clf:
    :param x:
    :param y:
    :param feature_names:
    :param nseq:
    :return:
    """
    # Get feature importances
    importances = clf.feature_importances_
    indices = np.argsort(importances)[::-1]

    if nseq is None:
        nseq = np.arange(1, len(importances), 10)

    # Retrain a classifier with only the N most important features
    accs = []
    for n in nseq:
        best_features = feature_names[indices[0:n]]
        x2 = x[:, indices[0:n]]
        x2_tr, x2_te, y2_tr, y2_te = sklearn.model_selection.train_test_split(
            x2, y, test_size=0.25)
        clf2 = sklearn.ensemble.RandomForestClassifier(n_estimators=100)
        clf2.fit(x2_tr, y2_tr)
        pred = clf2.predict(x2_te)
        mypd = pd.DataFrame(data={"Predicted": pred, "GroundTruth": y2_te})
        conf_matrix = pd.crosstab(
            index=mypd["Predicted"],
            columns=mypd["GroundTruth"])
        accs.append({
            "n": n,
            "classifier": (clf2, clf2.score(x2_te, y2_te), best_features),
            "confusion_matrix": conf_matrix})

    return accs


@DeprecationWarning
def optimize_clf(cell_line, save=True):
    """
    This function trains a classifier and then optimizes it to use as few
    features as necessary. Returns the classifier that requires as few
    features as possible while still delivering an accuracy within 1%-point
    of the full classifier

    :param cell_line:
    :param save: Boolean. If True, the function tries to load a preexisting
        classifier and also saves a newly trained classifier to disk.
    :return: Tuple (classifier, accuracy, features used)
    """

    clf_fn = os.path.join(BASEDIR, "OrganoidClassifier_%s.pkl" % cell_line)
    if save and os.path.isfile(clf_fn):
        with open(clf_fn, "rb") as f:
            return pickle.load(f)

    clf, acc, feature_names, x, y = train_classifier(
        cell_line=cell_line, return_data=True)

    reduced_classifiers = select_minimal_features_for_classifier(
        clf, x, y, feature_names, np.arange(1, 50, 1))

    acc_diff = np.abs(np.array(
        [s["classifier"][1] for s in reduced_classifiers]) - acc)
    first_good_index = np.where(acc_diff <= 0.01)[0][0]
    best_clf = reduced_classifiers[first_good_index]["classifier"]

    if save:
        with open(clf_fn, "wb") as f:
            pickle.dump(best_clf, f, pickle.HIGHEST_PROTOCOL)

    return best_clf


def apply_to_plate(plate):
    """
    Trains a classifier on pos/neg controls and then applies it to
    all wells of a plate.

    :param plate:
    :return:
    """

    if len(plate) != 14:
        raise ValueError("'plate' should be a 14-char string")

    clf = train_classifier(plate[0:7], save=True)

    # Load layout
    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(LAYOUTDIR, "%s.xlsx" % layout_id))

    well_dir = os.path.join(FEATUREDIR, plate, "wells_normalized")
    wells = sorted(
        [well for well in os.listdir(well_dir) if
         well.startswith(plate)])

    dead_percentage = []
    for well in wells:
        well_fn = os.path.join(well_dir, well)
        with h5py.File(well_fn, "r") as h5handle:
            features = np.transpose(
                h5handle["features_%s" % FEATURETYPE][()])
            feature_names = h5handle["feature_names_%s" % FEATURETYPE][()]

        # Reduce features to required set
        features = features[:, np.in1d(feature_names, clf[2])]

        # Remove bad organoids
        bad_organoids = np.sum(np.isfinite(features), axis=1)
        features = features[bad_organoids != 0, :]

        # Impute missing values
        col_medians = np.nanmedian(features, axis=0)
        col_medians = np.stack([col_medians] * features.shape[0])
        features[~np.isfinite(features)] = col_medians[~np.isfinite(features)]

        # Apply classifier to wells
        prediction = clf[0].predict(features)
        prediction_prob = clf[0].predict_proba(features)
        index_live = np.where(clf[0].classes_ == "NEG")[0][0]
        index_dead = np.where(clf[0].classes_ == "POS")[0][0]
        prediction_prob_live = np.nanmedian(
            prediction_prob[:, index_live][
                prediction_prob[:, index_live] > 0.5])
        prediction_prob_dead = np.nanmedian(
            prediction_prob[:, index_dead][
                prediction_prob[:, index_dead] > 0.5])

        well_id = "".join(well.split("_")[1:3])
        entry = layout.loc[
            layout["Well_ID_384"] == well_id,
            ("Product.Name", "concentration")]
        entry["Product.Name"] = entry["Product.Name"].astype(str)
        entry["concentration"] = entry["concentration"].astype(str)
        entry["Well.ID"] = well_id
        entry["Num.Objects"] = features.shape[0]
        entry["Percent.Live"] = np.mean(prediction == "NEG")
        entry["Percent.Dead"] = np.mean(prediction == "POS")
        entry["Median.Certainty.Live"] = prediction_prob_live
        entry["Median.Certainty.Dead"] = prediction_prob_dead
        dead_percentage.append(entry)
    dead_percentage = pd.concat(dead_percentage)

    return dead_percentage


if __name__ == "__main__":
    cmd_args = sys.argv
    if len(cmd_args) < 3:
        print(
            "Usage: %s [APPLY_TO|TRAIN] <PLATE|CELL_LINE>" % cmd_args[0])
        sys.exit()
    if cmd_args[1] == "TRAIN":
        train_classifier(cmd_args[2], return_data=False, save=True)
    elif cmd_args[1] == "APPLY_TO":
        myplate = cmd_args[2]
        results = apply_to_plate(myplate)
        out_fn = os.path.join(
            BASEDIR, myplate, "%s_organoid_classification.h5"
            % myplate)
        if not os.path.isdir(os.path.dirname(out_fn)):
            os.makedirs(os.path.dirname(out_fn))
        results.to_hdf(out_fn, myplate)
