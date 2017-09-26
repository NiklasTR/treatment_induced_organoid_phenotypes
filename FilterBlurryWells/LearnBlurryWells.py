# This script trains a random forest classifier to detect blurry wells and
# removes these from the feature set

from __future__ import division
import os
import h5py
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

if __name__ == "__main__":
    # Set up directories
    feature_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
    basedir = "/Users/jansauer/Thesis/Projects/PROMISE/FilterBlurryWells"

    # Use mouse plates
    organoid_type = "M"
    if organoid_type == "M":
        basedir = os.path.join(basedir, "mouse")
    elif organoid_type == "D":
        basedir = os.path.join(basedir, "human")
    else:
        print("'organoid_type' must be either 'M' or 'D'")
        quit()

    plates = [
        s for s in os.listdir(feature_dir)
        if s.startswith(organoid_type)]
    plates.sort()

    # Load inception features
    # Feature type can currently only be 'inception'
    feature_type = "inception"
    print("Loading '%s' features ..." % feature_type)
    features = []
    wells = []
    metadata = []
    for plate in plates:
        fn = os.path.join(feature_dir, feature_type, plate,
                          "%s_features_%s.h5" % (plate, feature_type))
        with h5py.File(fn, "r") as h5handle:
            feats = h5handle["features"][()]
            w = h5handle["well_names"][()]
            md = h5handle["metadata"][()]

        features.append(feats)
        wells += [plate + "_" + s for s in w]
        metadata.append(md)

    features = np.concatenate(features, axis=0)
    wells = np.array(wells)

    # Sanity check to ensure all metadata entries are identical
    md_deviations = 0
    for ii in range(len(metadata)-1):
        md_deviations += np.sum(metadata[ii] != metadata[ii+1])
    if md_deviations != 0:
        raise Exception("Metadata does not match in features")

    # Remove nan entries
    # (There aren't any anymore, but this is a good safeguard to have
    # regardless)
    nan_features = np.isnan(features)
    nan_rows = np.sum(nan_features, axis=1)
    sel_features = features[np.where(nan_rows == 0)]
    sel_wells = wells[np.where(nan_rows == 0)]

    # Load blurry well annotation
    anno_fn = os.path.join(
        basedir, "annotation", "blurry_wells_anno_%s.dat" % organoid_type)
    with open(anno_fn, "r") as f:
        blurry_wells = f.readlines()
    blurry_wells.sort()
    blurry_wells = np.array([s.strip() for s in blurry_wells])

    is_blurry = np.in1d(sel_wells, blurry_wells)
    negative_wells = sel_wells[is_blurry]
    negative_features = sel_features[is_blurry, :]

    # Select randomly non-blurry wells from the left half (col <= 12) of the
    # plates (except for M001W01P008L07!)
    n_wells = np.sum(is_blurry)
    sel_cols = np.array([int(s.split("_")[2]) for s in sel_wells])
    platenames = np.array([s.split("_")[0] for s in sel_wells])
    indices = (sel_cols <= 12) * (platenames != "M001W01P008L07")
    possible_wells = sel_wells[np.where(indices)]
    possible_features = sel_features[np.where(indices)]
    selection = np.random.choice(
        a=range(len(possible_wells)), size=n_wells, replace=False)
    positive_wells = possible_wells[selection]
    positive_features = possible_features[selection, :]

    # Train random forest
    print("Training random forest classifier ...")
    X = np.concatenate((negative_features, positive_features), axis=0)
    Y = np.repeat(("Blurry", "Sharp"), n_wells)
    X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size=0.25)

    rfc = RandomForestClassifier(n_estimators=500)
    rfc.fit(X_train, Y_train)

    rfc_score = rfc.score(X_val, Y_val)
    print("Random forest validation accuracy: %s" % rfc_score)
    with open(os.path.join(basedir, "rfc_accuracy.txt"), "w") as f:
        f.write(str(rfc_score) + "\n")

    # Run prediction on all wells
    print("Running trained classifier on all wells ...")
    rfc_predict = rfc.predict(sel_features)

    # Save files
    print("Saving outputs ...")
    all_blurry_wells = wells[rfc_predict == "Blurry"]
    all_sharp_wells = wells[rfc_predict == "Sharp"]

    new_blurry_wells = np.array([
        s for s in all_blurry_wells
        if s not in blurry_wells])

    blurry_out_file = os.path.join(basedir, "blurry_wells_predicted.txt")
    blurry_out_diff_file = os.path.join(
        basedir, "unannotated_blurry_wells_predicted.txt")
    sharp_out_file = os.path.join(basedir, "focused_wells_predicted.txt")
    with open(blurry_out_file, "w") as f:
        for well in all_blurry_wells:
            f.write(well + "\n")
    with open(sharp_out_file, "w") as f:
        for well in all_sharp_wells:
            f.write(well + "\n")
    with open(blurry_out_diff_file, "w") as f:
        f.write("# These are blurry wells that did not come from the\n")
        f.write("# manual annotation. They are useful for visualization\n")
        f.write("# of the classifier's output.\n")
        for well in new_blurry_wells:
            f.write(well + "\n")

    # The classifier trains rather quickly, but keep it for posterity anyway
    with open(os.path.join(basedir, "random_forest_classifier.pkl"), "w") as f:
        pickle.dump(rfc, f)

    # Save features of "good" wells in one place
    features_sharp = sel_features[np.where(rfc_predict == "Sharp")]
    features_blurry = sel_features[np.where(rfc_predict == "Blurry")]
    wells_sharp = sel_wells[np.where(rfc_predict == "Sharp")]
    wells_blurry = sel_wells[np.where(rfc_predict == "Blurry")]

    # Augment metadata
    metadata = metadata[0].tolist()
    metadata.append(
        ["Post-Processing",
         "Blurry wells removed with a random forest classifier"])

    features_out_fname = os.path.join(
        basedir, "Features_Cleaned_%s.h5" % organoid_type)

    with h5py.File(features_out_fname, "w-") as h5handle:
        h5handle.create_dataset(
            name="features_focused", data=features_sharp, compression=3)
        h5handle.create_dataset(
            name="features_blurry", data=features_blurry, compression=3)
        h5handle.create_dataset(
            name="well_names_focused", data=wells_sharp, compression=3)
        h5handle.create_dataset(
            name="well_names_blurry", data=wells_blurry, compression=3)
        h5handle.create_dataset(
            name="metadata", data=np.array(metadata), compression=3)
