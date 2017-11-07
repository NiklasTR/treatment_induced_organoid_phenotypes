import numpy as np
import Config
import os
import LoadFeatures
import NormalizeFeatures
import Utils
import h5py
import pandas as pd


def calc_feature_summary(
        features, feature_names,
        summary_func_middle,
        summary_func_var,
        kwargs_middle=None,
        kwargs_var=None,
        size_threshold=Config.SIZETHRESHOLD):
    """
    Calculate the summaries over features. Defaults to the median and MAD.
    This function splits the individual objects into organoids and shrapnel
    and treats them individually.

    :param features: A 2D numpy matrix with the shape (features, samples)
    :param feature_names:
    :param summary_func_middle:
    :param summary_func_var:
    :param kwargs_middle:
    :param kwargs_var:
    :param size_threshold:

    :return:
        - A feature vector for the well
        - A list of feature names
    """
    if kwargs_middle is None:
        kwargs_middle = dict()
    if kwargs_var is None:
        kwargs_var = {"center": np.nanmedian}

    # Remove the feature "FIELD"
    if "FIELD" in feature_names:
        features = np.delete(
            features, np.where(feature_names == "FIELD")[0], axis=0)
        feature_names = np.delete(
            feature_names, np.where(feature_names == "FIELD")[0])

    # Separate into organoids and shrapnel
    f_size = features[np.where(feature_names == "x.0.s.area")[0][0]]
    features_shrapnel = features[:, f_size < size_threshold]
    features_organoids = features[:, f_size >= size_threshold]

    # Shrapnel summaries
    if features_shrapnel.shape[1] == 0:
        features_shrapnel_m = np.zeros(
            shape=features_shrapnel.shape[0])
        features_shrapnel_v = np.zeros(
            shape=features_shrapnel.shape[0])
    else:
        features_shrapnel_m = summary_func_middle(
            features_shrapnel, axis=1, **kwargs_middle)
        features_shrapnel_v = summary_func_var(
            features_shrapnel, axis=1, **kwargs_var)
    num_shrapnel = features_shrapnel.shape[1]

    # Organoid summaries
    if features_organoids.shape[1] == 0:
        features_organoids_m = np.zeros(
            shape=features_organoids.shape[0])
        features_organoids_v = np.zeros(
            shape=features_organoids.shape[0])
    else:
        features_organoids_m = summary_func_middle(
            features_organoids, axis=1, **kwargs_middle)
        features_organoids_v = summary_func_var(
            features_organoids, axis=1, **kwargs_var)
    num_organoids = features_organoids.shape[1]

    return {
        "organoids": {
            "expected": features_organoids_m,
            "variation": features_organoids_v,
            "num_objects": num_organoids,
            "feature_names_expected": np.array([
                "organoids_%s_expected" % fn for fn in feature_names]),
            "feature_names_variation": np.array([
                "organoids_%s_variation" % fn for fn in feature_names])},
        "shrapnel": {
            "expected": features_shrapnel_m,
            "variation": features_shrapnel_v,
            "num_objects": num_shrapnel,
            "feature_names_expected": np.array([
                "shrapnel_%s_expected" % fn for fn in feature_names]),
            "feature_names_variation": np.array([
                "shrapnel_%s_variation" % fn for fn in feature_names])}}


def calc_well_summaries(plate):
    """
    Average the normalized organoid features on a plate separately for
    organoids and shrapnel

    :param plate:

    :return:
    """
    out_fn = os.path.join(
        Config.FEATUREDIR, plate,
        "%s_averaged_features_%s.h5" % (plate, Config.FEATURETYPE))
    if os.path.isfile(out_fn):
        return False

    wells = sorted([
        well for well in os.listdir(os.path.join(
            Config.FEATUREDIR, plate, "wells_normalized"))])

    features_combined_all = []
    features_combined_all_names = []
    for well in wells:
        well_id = "_".join(well.split("_")[0:3])
        features = LoadFeatures.load_organoid_features(
            wells=[well_id], normalized=True)

        features_organoids = features["features"][
            ..., features["object_type"] == "Organoid"]
        features_shrapnel = features["features"][
            ..., features["object_type"] == "Shrapnel"]

        features_organoids_summary = calc_feature_summary(
            features=features_organoids,
            feature_names=features["feature_names"],
            summary_func_middle=Utils.trim_func,
            summary_func_var=Utils.trim_func,
            kwargs_middle={"func": np.nanmean, "percent": 0.05},
            kwargs_var={"func": np.nanstd, "percent": 0.05},
            size_threshold=-np.inf)["organoids"]

        features_shrapnel_summary = calc_feature_summary(
            features=features_shrapnel,
            feature_names=features["feature_names"],
            summary_func_middle=Utils.trim_func,
            summary_func_var=Utils.trim_func,
            kwargs_middle={"func": np.nanmean, "percent": 0.05},
            kwargs_var={"func": np.nanstd, "percent": 0.05},
            size_threshold=np.inf)["shrapnel"]

        # Load DMSO averages for the number of objects
        dmso_summary = NormalizeFeatures.get_dmso_average_for_plate(plate)
        dmso_mean_num_organoids = dmso_summary["organoids"]["expected"][
            dmso_summary["organoids"]["feature_names_expected"] ==
            "organoids_num.of.objects_expected"]
        dmso_std_num_organoids = dmso_summary["organoids"]["variation"][
            dmso_summary["organoids"]["feature_names_variation"] ==
            "organoids_num.of.objects_variation"]
        dmso_mean_num_shrapnel = dmso_summary["shrapnel"]["expected"][
            dmso_summary["shrapnel"]["feature_names_expected"] ==
            "shrapnel_num.of.objects_expected"]
        dmso_std_num_shrapnel = dmso_summary["shrapnel"]["variation"][
            dmso_summary["shrapnel"]["feature_names_variation"] ==
            "shrapnel_num.of.objects_variation"]

        # Combine the expected and variation into a single array
        features_organoids_combined = np.concatenate((
            features_organoids_summary["expected"],
            features_organoids_summary["variation"],
            (features_organoids_summary["num_objects"] -
             dmso_mean_num_organoids) / dmso_std_num_organoids))
        features_shrapnel_combined = np.concatenate((
            features_shrapnel_summary["expected"],
            features_shrapnel_summary["variation"],
            (features_shrapnel_summary["num_objects"] -
             dmso_mean_num_shrapnel) / dmso_std_num_shrapnel))
        features_combined = np.concatenate((
            features_organoids_combined,
            features_shrapnel_combined))

        features_organoids_combined_names = np.concatenate((
            features_organoids_summary["feature_names_expected"],
            features_organoids_summary["feature_names_variation"],
            np.array(["organoids_num.of.objects"])))
        features_shrapnel_combined_names = np.concatenate((
            features_shrapnel_summary["feature_names_expected"],
            features_shrapnel_summary["feature_names_variation"],
            np.array(["shrapnel_num.of.objects"])))
        features_combined_names = np.concatenate((
            features_organoids_combined_names,
            features_shrapnel_combined_names))

        # Calculate Biomass in well
        seg_fn = os.path.join(
            Config.SEGMENTATIONDIR, plate,
            "%s_DNNsegmentation.h5" % well_id)
        with h5py.File(seg_fn, "r") as h5handle:
            mask = h5handle["mask"][()]
            well_biomass = np.sum(mask > 0)
        features_combined = np.append(features_combined, well_biomass)
        features_combined_names = np.append(
            features_combined_names, np.array(["Total.Biomass"]))

        features_combined_all.append(features_combined)
        features_combined_all_names.append(features_combined_names)

    features_combined_all = np.stack(features_combined_all, axis=1)
    fname_iter = iter(features_combined_all_names)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    features_combined_all_names = features_combined_all_names[0]

    wells = ["_".join(well.split("_")[0:3]) for well in wells]

    # Normalize the biomass by the DMSO summaries
    layout_id = plate[11:14]
    layout = pd.read_excel(
        io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))
    dmso_wells = layout.loc[
        layout["Product.Name"] == "DMSO",
        "Well_ID_384"].values
    dmso_wells = [
        "%s_%s_%s" % (plate, str(dmso_well[0]), str(dmso_well[1:3]))
        for dmso_well in dmso_wells]
    dmso_biomass = features_combined_all[
        features_combined_all_names == "Total.Biomass",
        np.in1d(wells, dmso_wells)]
    dmso_biomass_expected = Utils.trim_func(
        a=np.expand_dims(dmso_biomass, axis=0),
        func=np.nanmean, axis=1, percent=0.05)
    dmso_biomass_variation = Utils.trim_func(
        a=np.expand_dims(dmso_biomass, axis=0),
        func=np.nanstd, axis=1, percent=0.05)
    adj_biomass = (
        features_combined_all[features_combined_all_names == "Total.Biomass"] -
        dmso_biomass_expected) / dmso_biomass_variation
    features_combined_all[
        features_combined_all_names == "Total.Biomass",
        ...] = adj_biomass

    with h5py.File(out_fn, "w") as h5handle:
        h5handle.create_dataset(
            name="features_%s" % Config.FEATURETYPE,
            data=features_combined_all)
        h5handle.create_dataset(
            name="feature_names_%s" % Config.FEATURETYPE,
            data=features_combined_all_names)
        h5handle.create_dataset(
            name="well_names",
            data=wells)

    return True


def calc_cell_line_features(cell_line):
    """
    Normalize the features for a given cell line
    :param cell_line:
    :return:
    """

    cl_dir = os.path.join(Config.BASEDIR, "features")
    if not os.path.isdir(cl_dir):
        os.makedirs(cl_dir)

    cl_fn = os.path.join(
        cl_dir, "%s_averaged_features_%s.h5"
                % (cell_line, Config.FEATURETYPE))
    if os.path.isfile(cl_fn):
        print("Cell line features file already exists for '%s'" % cell_line)
        return False

    cl_plates = sorted([
        plate for plate in os.listdir(Config.FEATUREDIR)
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
    cl_plates = [
        cl_plates[i] for i in
        range(len(cl_plates)) if use_plate[i]]

    # Sort the plates
    cl_plates = sorted(cl_plates, key=lambda x: x[9:14])

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
    cl_concentrations = []
    for plate in cl_plates:
        feature_fn = os.path.join(
            Config.FEATUREDIR, plate,
            "%s_averaged_features_%s.h5" % (plate, Config.FEATURETYPE))
        with h5py.File(feature_fn, "r") as h5handle:
            features = h5handle["features_%s" % Config.FEATURETYPE][()]
            feature_names = h5handle["feature_names_%s" % Config.FEATURETYPE][()]
            well_names = h5handle["well_names"][()]
        cl_replicates.append([replicates[plate]] * len(well_names))
        cl_features.append(features)
        cl_feature_names.append(feature_names)
        cl_well_names.append(well_names)

        layout_id = plate[11:14]
        layout = pd.read_excel(
            io=os.path.join(Config.LAYOUTDIR, "%s.xlsx" % layout_id))
        well_ids = np.array([
            "".join(well_name.split("_")[1:3])
            for well_name in well_names])
        for well_id in well_ids:
            cl_drugs.append(layout.loc[
                                layout["Well_ID_384"] == well_id,
                                "Product.Name"].values)
            if "concentration" in layout.columns:
                cl_concentrations.append(layout.loc[
                                             layout["Well_ID_384"] == well_id,
                                             "concentration"].values)
            else:
                cl_concentrations.append(np.array([np.nan]))
    fname_iter = iter(cl_feature_names)
    if not all(np.array_equal(next(fname_iter), rest) for rest in fname_iter):
        raise Warning("Not all feature names were identical between files")
    cl_feature_names = cl_feature_names[0]
    cl_features = np.concatenate(cl_features, axis=1)
    cl_well_names = np.concatenate(cl_well_names, axis=0)
    cl_replicates = np.concatenate(cl_replicates, axis=0)
    cl_drugs = np.concatenate(cl_drugs, axis=0).astype(np.str)
    cl_concentrations = np.concatenate(cl_concentrations, axis=0)

    # Calculate z score across cell line and per replicate
    features_masked = np.ma.array(
        cl_features, mask=np.isnan(cl_features))
    # Z score across entire cell line
    features_z = scipy.stats.mstats.zscore(
        a=features_masked, axis=1)
    # Z score across replicate
    features_z_rep = np.zeros(
        shape=cl_features.shape)
    features_z_rep[..., cl_replicates == 1] = scipy.stats.mstats.zscore(
        a=features_masked[..., cl_replicates == 1], axis=1)
    features_z_rep[..., cl_replicates == 2] = scipy.stats.mstats.zscore(
        a=features_masked[..., cl_replicates == 2], axis=1)
    features_z_rep = features_z_rep

    # Save data
    with h5py.File(cl_fn, "w-") as h5handle:
        h5handle.create_dataset(
            name="features_zscore_%s" % Config.FEATURETYPE,
            data=features_z)
        h5handle.create_dataset(
            name="features_zscore_rep_%s" % Config.FEATURETYPE,
            data=features_z_rep)
        h5handle.create_dataset(
            name="features_%s" % Config.FEATURETYPE,
            data=cl_features)
        h5handle.create_dataset(
            name="feature_names",
            data=cl_feature_names)
        h5handle.create_dataset(name="well_names", data=cl_well_names)
        h5handle.create_dataset(name="drugs", data=cl_drugs)
        h5handle.create_dataset(name="replicates", data=cl_replicates)
        h5handle.create_dataset(
            name="concentrations", data=cl_concentrations)