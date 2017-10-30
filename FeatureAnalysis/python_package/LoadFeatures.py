import os
import Config
import h5py
import numpy as np
import re
import BlurryOrganoidClassifier


def load_organoid_features(wells=None, plates=None, normalized=False):
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
        "organoids": ("features_organoids", "feature_names_organoids"),
        "clumps": ("features_clumps", "feature_names_clumps")}

    if Config.FEATURETYPE not in keymap.keys():
        raise KeyError("'Config.FEATURETYPE' must be one of '%s'"
                       % str(keymap.keys()))

    hdf5_keys = keymap[Config.FEATURETYPE]

    # Set wells to load
    if wells is None and plates is None:
        raise ValueError("At least one of 'wells' and 'plates' must be a list")
    # If 'plates' is set but 'wells' is not
    elif wells is None and plates is not None:
        wells = sorted([
            s for plate in plates for s in
            os.listdir(os.path.join(Config.FEATUREDIR, plate, "wells"))])
    # If 'wells' is set but 'plates' is not
    elif wells is not None and plates is None:
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
    for well in wells:
        feature_fn = os.path.join(
            Config.FEATUREDIR, well.split("_")[0], well_folder, well)
        try:
            with h5py.File(feature_fn, "r") as h5handle:
                well_features = h5handle[hdf5_keys[0]][()]
                well_feature_names = h5handle[hdf5_keys[1]][()]
            well_names.append(np.repeat(
                re.sub("_[0-9a-zA-Z]*\.h5", "", well),
                well_features.shape[1]))
            features.append(well_features)
            feature_names.append(well_feature_names)
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

    return {
        "features": features, "feature_names": feature_names,
        "well_names": well_names}


def remove_blurry_organoids(features, feature_names, well_names):
    """
    Removes blurry organoids

    :param features: A 2D numpy array with the shape (features, samples)
    :param feature_names:
    :param well_names:
    :return:
    """

    clf = BlurryOrganoidClassifier.get_blurry_organoid_classifier(
        Config.BLURRYORGANOIDCLF)

    features_clf = np.array(features).transpose()
    features_clf = features_clf[
        ..., [f in clf["feature_names"] for f in feature_names]]

    is_focused = clf["clf"].predict(features_clf)
    features = features[..., [s == 1 for s in is_focused]]
    well_names = well_names[[s == 1 for s in is_focused]]

    return {
        "features": features, "feature_names": feature_names,
        "well_names": well_names}