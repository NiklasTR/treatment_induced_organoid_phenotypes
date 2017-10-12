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
import pandas as pd
import scipy.stats.mstats


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
    middle_func = np.nanmedian
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
        if "c" in kwargs.keys():
            c = kwargs["c"]
        else:
            c = 0.05
        return np.log((features + np.sqrt(features ** 2 + c ** 2)) / 2.0)
    else:
        raise ValueError("Unsupported transform_type ('%s')" % transform_type)


def calc_plate_median(plate, feature_dir, feature_type):
    """
    Calculate and return the plate summary
    :param plate:
    :param feature_dir:
    :param feature_type:
    :return:
    """
    wells = sorted([
        well for well in
        os.listdir(os.path.join(feature_dir, plate, "wells"))
        if well.startswith(plate) if well.endswith(".h5")])
    avg_features = []
    avg_feature_names = []
    for well in wells:
        well_fn = os.path.join(feature_dir, plate, "wells", well)
        try:
            with h5py.File(well_fn, "r") as h5handle:
                features = h5handle["features_%s" % feature_type][()]
                feature_names = h5handle["feature_names_%s" % feature_type][()]
            avg = calc_well_average(
                features=features, feature_names=feature_names,
                blurry_organoid_clf=blurry_organoid_clf)
            avg_features.append(avg[0])
            avg_feature_names.append(avg[1])
        except:
            avg_features.append(np.repeat(np.nan, 6286))

    fname_iter = iter(avg_feature_names)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    avg_feature_names = avg_feature_names[0]
    avg_features = np.stack(avg_features, axis=1)
    well_names = np.array([well[0:19] for well in wells])
    return {
        "features_median": avg_features,
        "feature_names": avg_feature_names,
        "well_names": well_names}


if __name__ == "__main__":
    basedir = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/preprocessing"
    feature_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
    layout_dir = "/collab-ag-fischer/PROMISE/layouts/python_friendly"
    segmentation_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/segmentation"
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

    all_plates = sorted([
        plate for plate in os.listdir(feature_dir) if
        plate.startswith("D0") or plate.startswith("M001")])

    # Go through each plate and:
    # - calculate the well summaries
    # - Apply glog transformation
    # - Subtract the DMSO median of each plate to correct for plate effects
    # The intermediate steps are stored separately to allow variations in the
    # normalization to be tested and the well median calculation is a
    # relatively lengthy process and unlikely to change, making it easier to
    # adapt to different normalization methods (e.g. replacing glog)
    for plate in all_plates:
        out_fn = os.path.join(
            feature_dir, plate,
            "%s_averaged_features_%s.h5" % (plate, feature_type))

        # Calc median summaries
        calc_summary = True
        if os.path.isfile(out_fn):
            with h5py.File(out_fn, "r") as h5handle:
                if "features_median" in h5handle.keys():
                    print("Medians for plate '%s' already processed" % plate)
                    calc_summary = False
        if calc_summary:
            print("Calculating medians for plate '%s' ..." % plate)
            median_dat = calc_plate_median(
                plate=plate, feature_dir=feature_dir,
                feature_type=feature_type)
            with h5py.File(out_fn, "w-") as h5handle:
                h5handle.create_dataset(
                    name="features_median", data=median_dat["features_median"])
                h5handle.create_dataset(
                    name="feature_names", data=median_dat["feature_names"])
                h5handle.create_dataset(
                    name="well_names", data=median_dat["well_names"])

        # Add the total biomass as a feature prior to transformation
        with h5py.File(out_fn, "r") as h5handle:
            avg_features = h5handle["features_median"][()]
            feature_names = h5handle["feature_names"][()]
            well_names = h5handle["well_names"][()]
        if "Total.Biomass" in feature_names:
            print("Biomass feature already calculated for plate '%s'" % plate)
        else:
            print("Calculating total biomass for plate '%s'" % plate)
            plate_biomass = []
            for well_name in well_names:
                seg_fn = os.path.join(
                    segmentation_dir, plate,
                    "%s_DNNsegmentation.h5" % well_name)
                with h5py.File(seg_fn, "r") as h5handle:
                    mask = h5handle["mask"][()]
                    plate_biomass.append(np.sum(mask > 0))
            avg_features = np.concatenate(
                (avg_features, np.expand_dims(plate_biomass, 0)))
            feature_names = np.concatenate(
                (feature_names, np.array(["Total.Biomass"])))
            with h5py.File(out_fn, "r+") as h5handle:
                del h5handle["features_median"]
                del h5handle["feature_names"]
                h5handle.create_dataset(
                    name="features_median", data=avg_features)
                h5handle.create_dataset(
                    name="feature_names", data=feature_names)

        # Transform
        with h5py.File(out_fn, "r+") as h5handle:
            if "features_glog" in h5handle.keys():
                print("glog transform for plate '%s' already calculated" % plate)
            else:
                print("Calculating glog transform for plate '%s' ..." % plate)
                avg_features = h5handle["features_median"][()]
                avg_features_glog = transform_features(
                    avg_features, "glog", c=0.05)
                h5handle.create_dataset(
                    name="features_glog", data=avg_features_glog)

        # Subtract the median DMSO controls
        calc_dmso_norm = True
        with h5py.File(out_fn, "r") as h5handle:
            if "features" in h5handle.keys():
                print("DMSO Median normalization for plate "
                      "'%s' already calculated" % plate)
                calc_dmso_norm = False

        if calc_dmso_norm:
            print("Calculate DMSO median normalization for plate '%s' ..." % plate)
            with h5py.File(out_fn, "r") as h5handle:
                avg_features = h5handle["features_glog"][()]
                well_names = h5handle["well_names"][()]
            well_ids = np.array([
                "".join(well_name.split("_")[1:3])
                for well_name in well_names])
            layout_id = plate[11:14]
            layout = pd.read_excel(
                io=os.path.join(layout_dir, "%s.xlsx" % layout_id))
            dmso_wells = layout.loc[
                layout["Product.Name"] == "DMSO",
                "Well_ID_384"].values
            dmso_features = avg_features[:, np.in1d(well_ids, dmso_wells)]
            median_dmso_features = np.nanmedian(dmso_features, axis=1)
            plate_norm_features = np.transpose(
                avg_features.transpose() - median_dmso_features)

            with h5py.File(out_fn, "r+") as h5handle:
                h5handle.create_dataset(
                    name="features",
                    data=plate_norm_features)
                h5handle.create_dataset(
                    name="readme", data=np.array(
                        ["Well Median -> glog -> Subtract DMSO median"]))

    # Go through each cell line and:
    # - Calculate the z score
    all_cell_lines = sorted(set([
        plate[0:7] for plate in os.listdir(feature_dir) if
        plate.startswith("D0") or plate.startswith("M001")]))

    for cell_line in all_cell_lines:
        cl_fn = os.path.join(
            basedir, "%s_averaged_features_%s.h5"
                     % (cell_line, feature_type))
        if os.path.isfile(cl_fn):
            print("Cell line features file already exists for '%s'" % cell_line)
            continue

        print("Combining cell line '%s'" % cell_line)

        cl_plates = sorted([
            plate for plate in all_plates
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
        cl_plates = [cl_plates[i] for i in range(len(cl_plates)) if use_plate[i]]

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
        for plate in cl_plates:
            feature_fn = os.path.join(
                feature_dir, plate,
                "%s_averaged_features_%s.h5" % (plate, feature_type))
            with h5py.File(feature_fn, "r") as h5handle:
                features = h5handle["features"][()]
                feature_names = h5handle["feature_names"][()]
                well_names = h5handle["well_names"][()]
            cl_replicates.append([replicates[plate]] * len(well_names))
            cl_features.append(features)
            cl_feature_names.append(feature_names)
            cl_well_names.append(well_names)

            layout_id = plate[11:14]
            layout = pd.read_excel(
                io=os.path.join(layout_dir, "%s.xlsx" % layout_id))
            well_ids = np.array([
                "".join(well_name.split("_")[1:3])
                for well_name in well_names])
            for well_id in well_ids:
                cl_drugs.append(layout.loc[
                    layout["Well_ID_384"] == well_id,
                    "Product.Name"].values)
        fname_iter = iter(cl_feature_names)
        if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
            raise Warning("Not all feature names were identical between files")
        cl_feature_names = cl_feature_names[0]
        cl_features = np.concatenate(cl_features, axis=1)
        cl_well_names = np.concatenate(cl_well_names, axis=0)
        cl_replicates = np.concatenate(cl_replicates, axis=0)
        cl_drugs = np.concatenate(cl_drugs, axis=0).astype(np.str)

        # Calculate z score
        cl_features_z = scipy.stats.mstats.zscore(a=cl_features, axis=1)

        # Save data
        with h5py.File(cl_fn, "w-") as h5handle:
            h5handle.create_dataset(name="features", data=cl_features_z)
            h5handle.create_dataset(
                name="feature_names", data=cl_feature_names)
            h5handle.create_dataset(name="well_names", data=cl_well_names)
            h5handle.create_dataset(name="drugs", data=cl_drugs)
            h5handle.create_dataset(name="replicates", data=cl_replicates)
