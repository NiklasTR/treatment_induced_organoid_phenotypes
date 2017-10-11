# The feature files are large, approximately 600 MB per plate
#   -> ca. 7GB per mouse cell line
#   -> ca. 3GB per human donor
# To make them more manageable, they are merged into a single file (per donor
# and cell line) and preprocessed. Preprocessing means:
#   - Rescaling the features.
#       Features are standardized via the robust median-z score:
#       z = (x - median(x)) / mad(x)
#   - Downsampling the data type. The raw features are in float64 format, but
#     float32 is completely sufficient in terms of accuracy

import h5py
import os
import numpy as np
import statsmodels.robust


def normalize_features(x):
    """
    Normalizes features to lie between -1 and 1 (max-min normalized).
    This method is preferable for several reasons:
      - The relation of the absolute values of features becomes irrelevant
      - Features are differently distributed, many of which are known to not be
        normally distributed (disqualifying mean standardization)
      - It is robust towards sparse features, e.g. features that are 0 in most
        but not all cases (disqualifying median standardization)
      - Neural networks usually work with max-min normalized features, making
        these features easier to compare later on.
    :param x:
    :return:
    """

    x_min = np.min(x, axis=1)
    x_max = np.max(x, axis=1)
    x_range = x_max - x_min
    # If x_range is 0 then the values are constant and x - min == 0, i.e.
    # it doesn't matter what it's divided by
    x_range[x_range == 0] = 1
    z = (x - np.expand_dims(x_min, 1)) / np.expand_dims(x_range, 1)
    return (2 * z) - 1


def get_averaged_features(plate_id, feature_type, feature_dir):
    """
    Load features for a plate and calculate the median and MAD values.

    Separates into large objects ("signal") and small objects ("noise")
    :param plate_id:
    :param feature_type:
    :param feature_dir:
    :return:
    """
    keymap = {
        "organoids": ("features", "feature_names", "well_names"),
        "clumps": ("features_clumps", "feature_names_clumps", "well_names_clumps"),
        "noseg": ("features_noseg", "feature_names_noseg", "well_names_noseg")}

    if feature_type not in keymap.keys():
        raise KeyError("'feature_type' must be one of '%s'"
                       % str(keymap.keys()))

    hdf5_keys = keymap[feature_type]

    hdf5_in_file = os.path.join(
        feature_dir, plate_id, "%s_features.h5" % plate_id)
    with h5py.File(hdf5_in_file, "r") as h5handle:
        features = h5handle[hdf5_keys[0]][()]
        feature_names = h5handle[hdf5_keys[1]][()]
        well_names = h5handle[hdf5_keys[2]][()]

    # Remove the feature "FIELD"
    features = np.delete(
        features, np.where(feature_names == "FIELD")[0], axis=0)

    # Calculate median and MAD and the number of all objects
    features_noise_summary = []
    features_signal_summary = []
    sorted_well_names = sorted(set(well_names))
    for well_name in sorted_well_names:
        well_features = features[:, well_names == well_name]

        # Separate objects by a threshold. This is arbitrarily set to be an
        # area of 5000 pixels. Anything smaller than this is considered
        # "noise" with potential use and anything larger than this is
        # "signal"
        size_threshold = 5000
        f_size = well_features[np.where(feature_names == "x.0.s.area")[0]][0]
        features_noise = well_features[:, f_size < size_threshold]
        features_signal = well_features[:, f_size >= size_threshold]

        features_noise_median = np.median(features_noise, axis=1)
        features_noise_mad = statsmodels.robust.mad(features_noise, axis=1, c=1)
        features_noise_summary.append(np.concatenate(
            (features_noise_median, features_noise_mad,
             (features_noise.shape[1],))))
        features_signal_median = np.median(features_signal, axis=1)
        features_signal_mad = statsmodels.robust.mad(features_signal, axis=1, c=1)
        features_signal_summary.append(np.concatenate(
            (features_signal_median, features_signal_mad,
             (features_signal.shape[1],))))

    features_noise_summary = np.stack(features_noise_summary, axis=1)
    features_signal_summary = np.stack(features_signal_summary, axis=1)

    new_feature_names = np.delete(
        feature_names, np.where(feature_names == "FIELD")[0])
    full_feature_names = np.concatenate((
        np.core.defchararray.add(new_feature_names, "_median"),
        np.core.defchararray.add(new_feature_names, "_mad"),
        np.array(["num.of.objects"])))

    return {
        "signal": features_signal_summary, "noise": features_noise_summary,
        "feature_names": full_feature_names, "well_names": sorted_well_names,
        "plate": plate_id}


def merge_all_averaged_features(cell_line, feature_type, feature_dir, out_dir):
    """

    :param cell_line:
    :param feature_type:
    :param feature_dir:
    :param out_dir:
    :return:
    """
    hdf5_out_fname = os.path.join(
        out_dir, "%s_features_averaged_%s.h5" % (cell_line, feature_type))
    if os.path.isfile(hdf5_out_fname):
        print("WARNING! File '%s' already exists.\n"
              "Because merging HDF5 files is a very time-consuming process,\n"
              "please make sure the existing file can be overwritten and\n"
              "delete it manually." % hdf5_out_fname)
        return False

    keymap = {
        "organoids": ("features", "feature_names", "well_names"),
        "clumps": ("features_clumps", "feature_names_clumps",
                   "well_names_clumps"),
        "noseg": ("features_noseg", "feature_names_noseg", "well_names_noseg")}

    if feature_type not in keymap.keys():
        raise KeyError("'feature_type' must be one of '%s'" % str(keymap.keys()))

    # Get all merged features
    all_plates = [
        plate for plate in os.listdir(feature_dir)
        if plate.startswith(cell_line)
        if os.path.isdir(os.path.join(feature_dir, plate))]

    features = []
    for plate in all_plates:
        features.append(get_averaged_features(
            plate_id=plate, feature_type=feature_type,
            feature_dir=feature_dir))

    # Test that feature names are identical
    all_feature_names = [tuple(f["feature_names"]) for f in features]
    if len(set(all_feature_names)) > 1:
        raise Exception("The features for the plates aren't identical")
    all_feature_names = all_feature_names[0]

    # Combine features into a single array
    # Well names
    all_well_names = [
        "%s_%s" % (f["plate"], w) for f in features
        for w in f["well_names"]]

    # Signal
    signal_features = np.concatenate([f["signal"] for f in features], axis=1)
    signal_features = signal_features.astype(np.float32)
    noise_features = np.concatenate([f["noise"] for f in features], axis=1)
    noise_features = noise_features.astype(np.float32)

    with h5py.File(hdf5_out_fname, "w-") as h5handle:
        h5handle.create_dataset(name="features_signal", data=signal_features)
        h5handle.create_dataset(name="features_noise", data=noise_features)
        h5handle.create_dataset(name="feature_names", data=all_feature_names)
        h5handle.create_dataset(name="well_names", data=all_well_names)
        h5handle.create_dataset(name="feature_type", data=feature_type)
        h5handle.create_dataset(name="cell_line", data=cell_line)

    return True


def load_hdf5_file_for_plate(plate_id, feature_type, feature_dir):
    """
    A helper function to load a feature type for a given plate
    :param plate_id:
    :param feature_type: Can be "organoids", "clumps", "noseg"
    :param feature_dir:
    :return: A Pandas data frame
    """
    fname = os.path.join(feature_dir, plate_id, "%s_features.h5" % plate_id)

    keymap = {
        "organoids": ("features", "feature_names", "well_names"),
        "clumps": ("features_clumps", "feature_names_clumps",
                   "well_names_clumps"),
        "noseg": ("features_noseg", "feature_names_noseg", "well_names_noseg")}

    if feature_type not in keymap.keys():
        raise KeyError("'feature_type' must be one of '%s'" % str(keymap.keys()))

    hdf5_keys = keymap[feature_type]

    with h5py.File(fname, "r") as h5handle:
        features = h5handle[hdf5_keys[0]][()]
        feature_names = h5handle[hdf5_keys[1]][()]
        well_names = h5handle[hdf5_keys[2]][()]

    # unique-ify the well names for easier recalling of indices
    # in the big file
    unique_wells_counts = dict(zip(*np.unique(
        well_names, return_counts=True)))
    well_ids = []
    for well in sorted(unique_wells_counts.keys()):
        well_ids += [
            "%s_%s" % (well, str(i)) for i
            in range(unique_wells_counts[well])]

    # Sanity check: Well names should stay unchanged
    if len(well_names) != len(well_ids):
        raise Exception("Critical Error @ well ID assignment")
    equal = [
        well_names[i] == well_ids[i][0:3]
        for i in range(len(well_names))]
    if not all(equal):
        raise Exception("Critical Error @ well ID assignment")

    # Add the plate ID to the well names for later identification
    well_names = np.core.defchararray.add("%s_" % plate_id, well_names)
    well_ids = np.core.defchararray.add("%s_" % plate_id, well_ids)

    return {
        "features": features, "feature_names": feature_names,
        "well_names": well_names, "well_ids": well_ids}


def create_full_feature_file(cell_line, feature_type, feature_dir, out_dir):
    """
    Merges all feature files for a given cell line into one hdf5 file.
    'cell_line' is the first 7 characters of the plate ID, e.g. 'M001A03' for
    mouse organoids or 'D004T01' for human donors
    :param cell_line:
    :param feature_type: Can be "organoids", "clumps", "noseg"
    :param feature_dir:
    :param out_dir:
    :return:
    """

    hdf5_out_fname = os.path.join(
        out_dir, "%s_features_%s.h5" % (cell_line, feature_type))
    if os.path.isfile(hdf5_out_fname):
        print("WARNING! File '%s' already exists.\n"
              "Because merging HDF5 files is a very time-consuming process,\n"
              "please make sure the existing file can be overwritten and\n"
              "delete it manually." % hdf5_out_fname)
        return False

    keymap = {
        "organoids": ("features", "feature_names", "well_names"),
        "clumps": ("features_clumps", "feature_names_clumps",
                   "well_names_clumps"),
        "noseg": ("features_noseg", "feature_names_noseg", "well_names_noseg")}

    if feature_type not in keymap.keys():
        raise KeyError("'feature_type' must be one of '%s'" % str(keymap.keys()))

    hdf5_keys = keymap[feature_type]

    # Go through all plates and determine the size of the output file
    all_plates = [
        plate for plate in os.listdir(feature_dir)
        if plate.startswith(cell_line)
        if os.path.isdir(os.path.join(feature_dir, plate))]

    feature_names = None
    well_names = []
    well_ids = []
    for plate in all_plates:
        hdf5_in_fname = os.path.join(
            feature_dir, plate, "%s_features.h5" % plate)
        with h5py.File(hdf5_in_fname, "r") as h5handle:
            plate_feature_names = h5handle[hdf5_keys[1]][()]
            plate_well_names = h5handle[hdf5_keys[2]][()]

        # unique-ify the well names for easier reference later
        unique_wells_counts = dict(zip(*np.unique(
            plate_well_names, return_counts=True)))
        plate_well_ids = []
        for well in sorted(unique_wells_counts.keys()):
            plate_well_ids += [
                "%s_%s" % (well, str(i)) for i
                in range(unique_wells_counts[well])]

        # Sanity check: Well names should stay unchanged
        if len(plate_well_names) != len(plate_well_ids):
            raise Exception("Critical Error @ well ID assignment")
        equal = [
            plate_well_names[i] == plate_well_ids[i][0:3]
            for i in range(len(plate_well_names))]
        if not all(equal):
            raise Exception("Critical Error @ well ID assignment")

        # Add the plate ID to the well names for later identification
        well_names.append(np.core.defchararray.add(
            "%s_" % plate, plate_well_names))
        well_ids.append(np.core.defchararray.add(
            "%s_" % plate, plate_well_ids))

        if feature_names is None:
            feature_names = plate_feature_names
        if not np.all(feature_names == plate_feature_names):
            raise Exception(
                "The features in '%s' don't match the expected "
                "features" % hdf5_in_fname)

    well_names = np.concatenate(well_names)
    well_ids = np.concatenate(well_ids)

    # Create the empty hdf5 file
    with h5py.File(hdf5_out_fname, "w-") as h5handle:
        h5handle.create_dataset(
            name="features", shape=(len(feature_names), len(well_names)),
            dtype=np.float32)
        h5handle.create_dataset(name="feature_names", data=feature_names)
        h5handle.create_dataset(name="well_names", data=well_names)
        h5handle.create_dataset(name="well_ids", data=well_ids)
        h5handle.create_dataset(name="feature_type", data=feature_type)

    # Go through each plate, like above, but this time load the features as well
    # and store them in the corresponding spots
    for plate in all_plates:
        plate_data = load_hdf5_file_for_plate(
            plate_id=plate, feature_type=feature_type,
            feature_dir=feature_dir)

        # Get the corresponding indices in the big file
        with h5py.File(hdf5_out_fname, "r") as h5handle:
            outfile_well_ids = h5handle["well_ids"][()]

        outfile_well_indices = np.where(np.in1d(
            outfile_well_ids, plate_data["well_ids"]))[0]

        # Slicing is faster than grabbing each individual index, so create
        # contiguous slices from the indices (this should be only one
        # continuous list)
        outfile_well_indices = np.array_split(
            outfile_well_indices,
            np.where(np.diff(outfile_well_indices) != 1)[0]+1)

        outfile_well_indices = [
            slice(min(entry), max(entry)+1) for entry in outfile_well_indices]

        with h5py.File(hdf5_out_fname, "r+") as h5handle:
            for entry in outfile_well_indices:
                h5handle["features"][:, entry] = plate_data["features"]

    return True


if __name__ == "__main__":
    featuredir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
    outdir = "/Users/jansauer/Thesis/Projects/PROMISE/FeatureAnalysis/features"
    all_cell_lines = sorted(set([
        p[0:7] for p in os.listdir(featuredir)
        if p.startswith("M")]))

    for cl in all_cell_lines:
        success = merge_all_averaged_features(
            cl, "organoids", featuredir, outdir)
        print cl, str(success)
