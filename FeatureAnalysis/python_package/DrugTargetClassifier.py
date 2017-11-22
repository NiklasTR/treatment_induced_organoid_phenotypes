# These functions try to classify organoids by the drug targets
from __future__ import division
import os
import Config
import pickle
import LoadFeatures
import pandas as pd
import numpy as np
import sklearn.model_selection
import sklearn.ensemble


def train_classifier(cell_line, return_data=False, save=True):
    """
    Train a classifier to organoids treated with different drugs. The
    drugs are grouped into their targets

    This function only uses layouts L02 and L03, i.e. without drug
    concentration data.

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

    clf_fn = os.path.join(
        Config.DRUGTARGETCLASSIFIERDIR,
        "OrganoidClassifier_%s.pkl" % cell_line)
    if save and os.path.isfile(clf_fn):
        with open(clf_fn, "rb") as f:
            return pickle.load(f)

    if not os.path.isdir(Config.DRUGTARGETCLASSIFIERDIR):
        os.makedirs(Config.DRUGTARGETCLASSIFIERDIR)

    # Set plates
    plates = sorted([
        plate for plate in os.listdir(Config.FEATUREDIR)
        if plate.startswith(cell_line) if plate.endswith("L02") or
        plate.endswith("L03")])

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

    # Load features
    cl_features = LoadFeatures.load_organoid_features(
        plates=plates, normalized=True)

    # Remove features with any non-finite values
    valid_features = np.sum(~np.isfinite(cl_features["features"]), axis=1) == 0
    cl_features["features"] = cl_features["features"][valid_features, ...]
    cl_features["feature_names"] = cl_features["feature_names"][valid_features]

    # Load the targets for each well
    cl_drugs = []
    for plate in plates:
        layout_id = plate[11:14]
        layout = pd.read_excel(
            io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))
        pl_wn = cl_features["well_names"][
            [well.startswith(plate) for well in cl_features["well_names"]]]
        pl_wn = [well.split("_")[1] + well.split("_")[2] for well in pl_wn]
        drugs = [
            layout.loc[layout["Well_ID_384"] == well,
                       "Product.Name"].values.astype(str)[0]
            for well in pl_wn]
        cl_drugs.append(drugs)
    cl_drugs = np.concatenate(cl_drugs)

    # Subsample the dataset to account for different numbers of organoids for
    # each treatment
    num_samples = np.unique(cl_drugs, return_counts=True)[1].min()
    x = []
    y = []
    for drug in np.unique(cl_drugs):
        drug_features = np.transpose(
            cl_features["features"][..., cl_drugs == drug])
        x.append(drug_features[np.random.choice(
            range(drug_features.shape[0]), num_samples, replace=False), ...])
        y.append(np.repeat(drug, num_samples))
    x = np.concatenate(x, axis=0)
    y = np.concatenate(y)
    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(
        x, y, test_size=0.25)

    # Run classifier
    clf = sklearn.ensemble.RandomForestClassifier(n_estimators=100)
    clf.fit(x_train, y_train)
    acc = clf.score(x_test, y_test)

    if save:
        with open(clf_fn, "wb") as f:
            pickle.dump((
                clf, acc, cl_features["feature_names"]),
                f, pickle.HIGHEST_PROTOCOL)

    if return_data:
        return clf, acc, cl_features["feature_names"], x, y
    else:
        return clf, acc, cl_features["feature_names"]
