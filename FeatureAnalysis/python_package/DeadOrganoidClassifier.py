import h5py
import Config
import numpy as np
import pandas as pd
import os
import pickle
import sklearn.ensemble
import sklearn.model_selection
import sklearn.metrics
import sklearn.preprocessing


def get_classification_labels(
        features, feature_names, well_names, object_type, return_probs=False):
    """
    Gets the dead/live label of a given set of features. Returns a (numpy)
    array of labels in which any original
    :param features:
    :param feature_names:
    :param well_names:
    :param object_type:
    :param return_probs:
    :return:
    """
    # Handle each cell line separately
    cell_lines = np.array([wn[0:7] for wn in well_names])
    out_vec = []
    for cl in set(cell_lines):
        cl_features = np.transpose(features[..., cell_lines == cl])

        # Impute missing features
        # This is done by calculating the mean and standard deviation of each
        # feature and assigning all NA entries randomly selected values from
        # a normal distribution. If all features are NA, then a distribution
        # N(0, 1) is used
        # The imputation is done for each training set individually
        masked_features = np.ma.array(
            data=cl_features,
            mask=~np.isfinite(cl_features))
        for ii in range(len(feature_names)):
            na_features = ~np.isfinite(cl_features[..., ii])
            if np.sum(na_features) == 0:
                continue
            if np.sum(~na_features) == 0:
                f_mean_pos = 0
                f_sd_pos = 1
            else:
                f_mean_pos = np.mean(masked_features[..., ii])
                f_sd_pos = np.std(masked_features[..., ii])
            cl_features[..., ii][na_features] = np.random.normal(
                f_mean_pos, f_sd_pos, np.sum(na_features))

        clf = train_classifier(cl, save=False)
        if return_probs:
            out_vec.append(clf[0].predict_proba(cl_features))
        else:
            out_vec.append(clf[0].predict(cl_features))

    out_vec = np.concatenate(out_vec, axis=0)
    return out_vec


def classify_organoids(plate):
    """
    Classifies the organoids into either live or dead organoids
    :param plate:
    :return:
    """

    if len(plate) != 14:
        raise ValueError("'plate' should be a 14-char string")

    out_fn = os.path.join(
        Config.DEADORGANOIDCLASSIFIERDIR,
        "%s_organoid_classification.csv" % plate)
    if os.path.isfile(out_fn):
        return pd.read_csv(out_fn)

    clf = train_classifier(plate[0:7], save=True)

    # Load layout
    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))

    well_dir = os.path.join(Config.FEATUREDIR, plate, "wells_normalized")
    wells = sorted(
        [well for well in os.listdir(well_dir) if
         well.startswith(plate)])

    dead_percentage = []
    for well in wells:
        well_fn = os.path.join(well_dir, well)
        with h5py.File(well_fn, "r") as h5handle:
            features = np.transpose(
                h5handle["features_%s" % Config.FEATURETYPE][()])
            feature_names = h5handle["feature_names_%s" % Config.FEATURETYPE][()]

        # Remove bad organoids
        bad_organoids = np.sum(np.isfinite(features), axis=1)
        features = features[bad_organoids != 0, :]

        # Impute missing features
        # This is done by calculating the mean and standard deviation of each
        # feature and assigning all NA entries randomly selected values from
        # a normal distribution. If all features are NA, then a distribution
        # N(0, 1) is used
        # The imputation is done for each training set individually
        masked_features = np.ma.array(
            data=features,
            mask=~np.isfinite(features))
        for ii in range(len(feature_names)):
            na_features = ~np.isfinite(features[..., ii])
            if np.sum(na_features) == 0:
                continue
            if np.sum(~na_features) == 0:
                f_mean_pos = 0
                f_sd_pos = 1
            else:
                f_mean_pos = np.mean(masked_features[..., ii])
                f_sd_pos = np.std(masked_features[..., ii])
            features[..., ii][na_features] = np.random.normal(
                f_mean_pos, f_sd_pos, np.sum(na_features))

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
    dead_percentage.to_csv(out_fn, index=False)

    return dead_percentage


def train_classifier(cell_line, save=True):
    """
    Train a classifier to learn the difference between a live and a dead
    organoid. The live and dead organoids are taken from "control" wells.

    Positive controls are, as determined visually and via initial feature
    analysis, Bortezomib and Irinotecan / SN-38 at concentrations of 0.2
    and 1.0. Negative controls are DMSO wells.

    :param cell_line:
    :param save: Boolean. If True, the function tries to load a preexisting
        classifier and also saves a newly trained classifier to disk.
    :return: A tuple (classifier, accuracy, features used) or
        (classifier, accuracy, features used, data_x, data_y) if
        return_data == True
    """
    if len(cell_line) != 7:
        raise ValueError("'cell_line' should be a 7-char string")

    clf_fn = os.path.join(
        Config.DEADORGANOIDCLASSIFIERDIR,
        "OrganoidClassifier_%s.pkl" % cell_line)
    if os.path.isfile(clf_fn):
        with open(clf_fn, "rb") as f:
            return pickle.load(f)

    if not os.path.isdir(Config.DEADORGANOIDCLASSIFIERDIR):
        os.makedirs(Config.DEADORGANOIDCLASSIFIERDIR)

    # Define controls
    pos_ctrls = [
        ("Bortezomib", (0.2, 1.0)),
        ("Irinotecan / SN-38", (0.2, 1.0))]
    neg_ctrls = [("DMSO", None)]

    # Run over one cell line
    plates = sorted([
        plate for plate in os.listdir(Config.FEATUREDIR)
        if plate.startswith(cell_line) if plate.endswith("L08")])

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
    plates = [
        plates[i] for i in
        range(len(plates)) if use_plate[i]]

    # Sort the plates
    plates = sorted(plates, key=lambda x: x[9:14])

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
                Config.FEATUREDIR, plate, "wells_normalized",
                "%s_%s_%s_features_normalized.h5" %
                (plate, well[0], well[1:3]))
            with h5py.File(features_fn, "r") as h5handle:
                pos_ctrl_features.append(
                    h5handle["features_%s" % Config.FEATURETYPE][()])
                pos_ctrl_feature_names.append(
                    h5handle["feature_names_%s" % Config.FEATURETYPE][()])
                pos_ctrl_object_types.append(
                    h5handle["object_type_%s" % Config.FEATURETYPE][()])

        for well in neg_ctrl_wells:
            features_fn = os.path.join(
                Config.FEATUREDIR, plate, "wells_normalized",
                "%s_%s_%s_features_normalized.h5" %
                (plate, well[0], well[1:3]))
            with h5py.File(features_fn, "r") as h5handle:
                neg_ctrl_features.append(
                    h5handle["features_%s" % Config.FEATURETYPE][()])
                neg_ctrl_feature_names.append(
                    h5handle["feature_names_%s" % Config.FEATURETYPE][()])
                neg_ctrl_object_types.append(
                    h5handle["object_type_%s" % Config.FEATURETYPE][()])

    # Check feature names
    pos_iter = iter(pos_ctrl_feature_names + neg_ctrl_feature_names)
    if not all(np.array_equal(next(pos_iter), rest) for rest in pos_iter):
        raise Warning("Not all feature names were identical between files")
    feature_names = pos_ctrl_feature_names[0]

    pos_ctrl_features = np.transpose(
        np.concatenate(pos_ctrl_features, axis=1))
    neg_ctrl_features = np.transpose(
        np.concatenate(neg_ctrl_features, axis=1))

    pos_ctrl_object_types = np.concatenate(pos_ctrl_object_types)
    neg_ctrl_object_types = np.concatenate(neg_ctrl_object_types)

    # Remove blurry organoids from training
    pos_ctrl_features = pos_ctrl_features[pos_ctrl_object_types != "BLURRY"]
    pos_ctrl_object_types = pos_ctrl_object_types[
        pos_ctrl_object_types != "BLURRY"]
    neg_ctrl_features = neg_ctrl_features[neg_ctrl_object_types != "BLURRY"]
    neg_ctrl_object_types = neg_ctrl_object_types[
        neg_ctrl_object_types != "BLURRY"]

    # Remove bad organoids from dataset (i.e. no features at all)
    bad_wells_pos = np.sum(np.isfinite(pos_ctrl_features), axis=1)
    pos_ctrl_features = pos_ctrl_features[bad_wells_pos >= 0, :]
    pos_ctrl_object_types = pos_ctrl_object_types[bad_wells_pos >= 0]
    bad_wells_neg = np.sum(np.isfinite(neg_ctrl_features), axis=1)
    neg_ctrl_features = neg_ctrl_features[bad_wells_neg >= 0, :]
    neg_ctrl_object_types = neg_ctrl_object_types[bad_wells_neg >= 0]

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

    # Impute missing features
    # This is done by assigning all missing values a random number from a
    # N(0, 1) distribution
    masked_x = np.ma.array(
        data=x,
        mask=~np.isfinite(x))
    x[masked_x.mask] = np.random.normal(
        0, 1, np.sum(masked_x.mask))

    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(
        x, y, test_size=0.25)

    # Run classifier
    clf = sklearn.ensemble.RandomForestClassifier(n_estimators=100)
    clf.fit(x_train, y_train)
    acc = clf.score(x_test, y_test)

    if save:
        with open(clf_fn, "wb") as f:
            pickle.dump((
                clf, acc, feature_names, x_test, y_test),
                f, pickle.HIGHEST_PROTOCOL)

    return clf, acc, feature_names, x_test, y_test


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
        io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))
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


def create_transfer_learning_data():
    """
    This function generates diagnostic data for the dead/live classifiers.
    This includes running each classifier on every other cell line's
    validation data.
    :return:
    """

    all_cell_lines = sorted(set([
        plate[0:7] for plate in os.listdir(Config.FEATUREDIR)
        if plate.startswith("D0")]))

    # Load validation data and classifiers
    val_data = []
    clfs = []
    for cl in all_cell_lines:
        clf = train_classifier(cl, save=False)
        val_data.append((clf[3], clf[4]))
        clfs.append((clf[0], clf[2]))

    # Calculate accuracy matrix
    acc_matrix = []
    for clf in clfs:
        entry = []
        for ii in range(len(val_data)):
            vd = val_data[ii]
            entry.append(clf[0].score(*vd))
        acc_matrix.append(entry)

    out_fn = os.path.join(
        Config.DEADORGANOIDCLASSIFIERDIR,
        "transfer_learning_matrix.csv")
    with open(out_fn, "w") as f:
        f.write("# Accuracy of classifier for cell line (ROW) applied to "
                "validation data for cell line (COL)\n")
        f.write("CLASSIFIER_DATA," + ",".join(all_cell_lines) + "\n")
        for ii in range(len(acc_matrix)):
            f.write(all_cell_lines[ii] + "," +
                    ",".join([str(s) for s in acc_matrix[ii]]) + "\n")


def create_roc_data(cell_line, data_cell_line=None, plot=False):
    """
    Create ROC data for the given cell line. Can optionally create ROC data
    using a classifier for one cell line and the validation data of another
    cell line.
    :param cell_line:
    :param data_cell_line:
    :param plot: Boolean. Determines whether data is returned or a plot shown

    :return:
    """
    if data_cell_line is None:
        data_cell_line = cell_line
    clf = train_classifier(cell_line, save=False)
    classifier = clf[0]
    val_clf = train_classifier(data_cell_line, save=False)
    x_val = val_clf[3]
    y_val = val_clf[4]

    y_true = sklearn.preprocessing.label_binarize(
        y_val, classes=classifier.classes_)[..., 0]
    y_pred = classifier.predict_proba(x_val)[..., 1]
    fpr, tpr, thresh = sklearn.metrics.roc_curve(y_true, y_pred)
    roc_auc = sklearn.metrics.auc(fpr, tpr)

    if plot:
        import matplotlib.pyplot as plt
        plt.figure()
        lw = 2
        plt.plot(fpr, tpr, color='darkorange',
                 lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic example')
        plt.legend(loc="lower right")
        plt.show()
    else:
        out_fn = os.path.join(
            Config.DEADORGANOIDCLASSIFIERDIR,
            "roc_data_%s_on_%s.csv" % (cell_line, data_cell_line))
        out_dat = np.stack((thresh, fpr, tpr), axis=1)
        with open(out_fn, "w") as f:
            f.write(
                "# ROC Data for %s classifier applied to %s data\n"
                % (cell_line, data_cell_line))
            f.write("# AUC: %s\n" % roc_auc)
            f.write("Threshold,FalsePosRate,TruePosRate\n")
            for ii in range(out_dat.shape[0]):
                f.write(",".join(out_dat[ii, ...].astype(str)) + "\n")


def train_classifier_on_all_cell_lines(save=True):
    """
    This version of the function trains a classifier on all cell lines rather
    than an individual one. This has the potential advantage of being able
    to account for all "dead" and "live" phenotypes that the controls might
    not cover.
    :param save:
    :return:
    """
    clf_fn = os.path.join(
        Config.DEADORGANOIDCLASSIFIERDIR, "full_classifier",
        "OrganoidClassifier_all_celllines.pkl")
    if os.path.isfile(clf_fn):
        with open(clf_fn, "rb") as f:
            return pickle.load(f)

    if not os.path.isdir(Config.DEADORGANOIDCLASSIFIERDIR):
        os.makedirs(Config.DEADORGANOIDCLASSIFIERDIR)

    # Define controls
    pos_ctrls = [
        ("Bortezomib", (0.2, 1.0)),
        ("Irinotecan / SN-38", (0.2, 1.0))]
    neg_ctrls = [("DMSO", None)]

    # Run over all cell lines
    plates = sorted([
        plate for plate in os.listdir(Config.FEATUREDIR)
        if plate.endswith("L08")])

    # Use only re-imaged plates (9XX vs 0XX)
    plate_ids = [s[0:4] + s[8:11] for s in plates]
    use_plate = []
    for plate_id in plate_ids:
        if plate_id[4] == "9":
            use_plate.append(True)
            continue
        reimaged = plate_id[0:4] + "9" + plate_id[5:7]
        if reimaged in plate_ids:
            use_plate.append(False)
        else:
            use_plate.append(True)
    plates = [
        plates[i] for i in
        range(len(plates)) if use_plate[i]]

    # Sort the plates
    plates = sorted(plates, key=lambda x: x[9:14])

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
                Config.FEATUREDIR, plate, "wells_normalized",
                "%s_%s_%s_features_normalized.h5" %
                (plate, well[0], well[1:3]))
            with h5py.File(features_fn, "r") as h5handle:
                pos_ctrl_features.append(
                    h5handle["features_%s" % Config.FEATURETYPE][()])
                pos_ctrl_feature_names.append(
                    h5handle["feature_names_%s" % Config.FEATURETYPE][()])
                pos_ctrl_object_types.append(
                    h5handle["object_type_%s" % Config.FEATURETYPE][()])

        for well in neg_ctrl_wells:
            features_fn = os.path.join(
                Config.FEATUREDIR, plate, "wells_normalized",
                "%s_%s_%s_features_normalized.h5" %
                (plate, well[0], well[1:3]))
            with h5py.File(features_fn, "r") as h5handle:
                neg_ctrl_features.append(
                    h5handle["features_%s" % Config.FEATURETYPE][()])
                neg_ctrl_feature_names.append(
                    h5handle["feature_names_%s" % Config.FEATURETYPE][()])
                neg_ctrl_object_types.append(
                    h5handle["object_type_%s" % Config.FEATURETYPE][()])

    # Check feature names
    pos_iter = iter(pos_ctrl_feature_names + neg_ctrl_feature_names)
    if not all(np.array_equal(next(pos_iter), rest) for rest in pos_iter):
        raise Warning("Not all feature names were identical between files")
    feature_names = pos_ctrl_feature_names[0]

    pos_ctrl_features = np.transpose(
        np.concatenate(pos_ctrl_features, axis=1))
    neg_ctrl_features = np.transpose(
        np.concatenate(neg_ctrl_features, axis=1))

    pos_ctrl_object_types = np.concatenate(pos_ctrl_object_types)
    neg_ctrl_object_types = np.concatenate(neg_ctrl_object_types)

    # Remove blurry organoids from training
    pos_ctrl_features = pos_ctrl_features[pos_ctrl_object_types != "BLURRY"]
    pos_ctrl_object_types = pos_ctrl_object_types[
        pos_ctrl_object_types != "BLURRY"]
    neg_ctrl_features = neg_ctrl_features[neg_ctrl_object_types != "BLURRY"]
    neg_ctrl_object_types = neg_ctrl_object_types[
        neg_ctrl_object_types != "BLURRY"]

    # Remove bad organoids from dataset (i.e. no features at all)
    bad_wells_pos = np.sum(np.isfinite(pos_ctrl_features), axis=1)
    pos_ctrl_features = pos_ctrl_features[bad_wells_pos >= 0, :]
    pos_ctrl_object_types = pos_ctrl_object_types[bad_wells_pos >= 0]
    bad_wells_neg = np.sum(np.isfinite(neg_ctrl_features), axis=1)
    neg_ctrl_features = neg_ctrl_features[bad_wells_neg >= 0, :]
    neg_ctrl_object_types = neg_ctrl_object_types[bad_wells_neg >= 0]

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

    # Impute missing features
    # Missing features are given a random number selected from a N(0, 1)
    # distribution
    masked_x = np.ma.array(
        data=x, mask=~np.isfinite(x))
    x[masked_x.mask] = np.random.normal(
            0, 1, np.sum(masked_x.mask))

    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(
        x, y, test_size=0.25)

    # Run classifier
    clf = sklearn.ensemble.RandomForestClassifier(n_estimators=100)
    clf.fit(x_train, y_train)
    acc = clf.score(x_test, y_test)

    if save:
        with open(clf_fn, "wb") as f:
            pickle.dump((
                clf, acc, feature_names, x_test, y_test),
                f, pickle.HIGHEST_PROTOCOL)

    return clf, acc, feature_names, x_test, y_test


def classify_organoids_with_full_classifier(plate):
    """
    Classifies the organoids into either live or dead organoids
    :param plate:
    :return:
    """

    if len(plate) != 14:
        raise ValueError("'plate' should be a 14-char string")

    out_fn = os.path.join(
        Config.DEADORGANOIDCLASSIFIERDIR, "full_classifier",
        "%s_organoid_classification_full_classifier.csv" % plate)
    if os.path.isfile(out_fn):
        return pd.read_csv(out_fn)

    clf = train_classifier_on_all_cell_lines(save=True)

    # Load layout
    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))

    well_dir = os.path.join(Config.FEATUREDIR, plate, "wells_normalized")
    wells = sorted(
        [well for well in os.listdir(well_dir) if
         well.startswith(plate)])

    dead_percentage = []
    for well in wells:
        well_fn = os.path.join(well_dir, well)
        with h5py.File(well_fn, "r") as h5handle:
            features = np.transpose(
                h5handle["features_%s" % Config.FEATURETYPE][()])
            feature_names = h5handle["feature_names_%s" % Config.FEATURETYPE][()]

        # Remove bad organoids
        bad_organoids = np.sum(np.isfinite(features), axis=1)
        features = features[bad_organoids != 0, :]

        # Impute missing features
        # Missing features are replaced with a random number from a N(0, 1)
        # distribution
        masked_features = np.ma.array(
            data=features,
            mask=~np.isfinite(features))
        features[masked_features.mask] = np.random.normal(
            0, 1, np.sum(masked_features.mask))

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
    dead_percentage.to_csv(out_fn, index=False)

    return dead_percentage


def create_transfer_learning_data_for_full_classifier():
    """
    This function generates diagnostic data for the dead/live classifiers.
    This includes running each classifier on every other cell line's
    validation data.
    :return:
    """

    all_cell_lines = sorted(set([
        plate[0:7] for plate in os.listdir(Config.FEATUREDIR)
        if plate.startswith("D0")]))

    # Load validation data and classifiers
    val_data = []
    clfs = [train_classifier_on_all_cell_lines(save=False)]
    for cl in all_cell_lines:
        clf = train_classifier(cl, save=False)
        val_data.append((clf[3], clf[4]))

    # Calculate accuracy matrix
    acc_matrix = [[]]
    for ii in range(len(val_data)):
        vd = val_data[ii]
        acc_matrix[0].append(clfs[0][0].score(*vd))

    out_fn = os.path.join(
        Config.DEADORGANOIDCLASSIFIERDIR, "full_classifier",
        "transfer_learning_matrix_full_classifier.csv")
    with open(out_fn, "w") as f:
        f.write("# Accuracy of full classifier (ROW) applied to "
                "validation data for cell line (COL)\n")
        f.write("CLASSIFIER_DATA," + ",".join(all_cell_lines) + "\n")
        for ii in range(len(acc_matrix)):
            f.write("Full Classifier" + "," +
                    ",".join([str(s) for s in acc_matrix[ii]]) + "\n")
