from __future__ import division

import random
import numpy as np
import pandas as pd
import h5py
import os
import scipy.ndimage
import sys
sys.path.append("/Users/jansauer/Thesis/Projects/PROMISE/DeepLearningFeatures")
import ImagePreprocessor


def get_wells_from_layout(layout_file, treatments, concentrations=None):
    """
    Get the IDs of plates treated with a given drug
    :param layout_file:
    :param treatments:
    :param concentrations:
    :return:
    """

    layout = pd.read_excel(layout_file)
    if concentrations is None:
        return layout.loc[np.in1d(
            layout["Product.Name"].values.astype(str),
            np.array(treatments)), "Well_ID_384"].values.astype(str)
    else:
        return layout.loc[np.logical_and(
            np.in1d(layout["Product.Name"].values.astype(str),
                    np.array(treatments)),
            np.in1d(layout["concentration"].values.astype(float),
                    concentrations)), "Well_ID_384"].values.astype(str)


def create_patches(
        plate, row, col, patch_rad, num_patches, promise_dir):
    """
    Load a well and select patches from the segmented foreground.
    Patches are guaranteed to contain at least 50% foreground in them.
    :param plate:
    :param row:
    :param col:
    :param patch_rad:
    :param num_patches:
    :param promise_dir:
    :return:
    """

    # Load images
    segmask_fn = os.path.join(
        promise_dir, "segmentation", plate,
        "%s_%s_%s_DNNsegmentation.h5" % (plate, row, col))
    with h5py.File(segmask_fn, "r") as h5handle:
        segmask = h5handle["mask"][()]
    segmask[segmask > 0] = 1

    proj_fn = os.path.join(
        promise_dir, "hdf5projection", plate,
        "%s_%s_%s_contrastProjections.h5" % (plate, row, col))
    with h5py.File(proj_fn, "r") as h5handle:
        proj = h5handle["images"][()]
    proj = np.stack([ImagePreprocessor.preprocess_image(img) for img in proj])

    # Apply mean filter to segmask and keep the pixels that correspond to
    # regions with >= 50% foreground
    segmask = np.stack([scipy.ndimage.uniform_filter(
        img.astype(np.float_), size=patch_rad) >= 0.5 for img in segmask])

    # Extract foreground points from reduced segmask
    points = np.transpose(np.where(segmask))

    # Remove points too close to the image borders
    valid_points_ax1 = np.logical_and(
        patch_rad <= points[:, 1],
        points[:, 1] < proj.shape[1] - patch_rad)
    valid_points_ax2 = np.logical_and(
        patch_rad <= points[:, 2],
        points[:, 2] < proj.shape[2] - patch_rad)
    valid_points = np.logical_and(valid_points_ax1, valid_points_ax2)
    points = points[valid_points]

    # Randomly sample points
    points_indices = np.random.choice(len(points), num_patches)
    points = points[points_indices]

    patches = np.zeros(shape=(len(points), patch_rad*2+1, patch_rad*2+1, 3))
    for point_index in range(len(points)):
        pp = points[point_index]
        patches[point_index, ...] = proj[
            pp[0], (pp[1] - patch_rad):(pp[1] + patch_rad + 1),
            (pp[2] - patch_rad):(pp[2] + patch_rad + 1), :]
    return patches


def create_training_data_pos_neg(
        patches_per_well=10000,
        patch_rad=40,
        promise_dir="/Users/jansauer/tmp/PROMISE_local",
        lib_dir="/Users/jansauer/tmp/watchdog_local/layouts",
        out_dir="."):
    """
    Create a training data hdf5 file from the libraries. This currently only
    uses the L08 plates, as only this library has positive and negative
    controls in known concentrations.
    :param patches_per_well:
    :param patch_rad:
    :param promise_dir:
    :param lib_dir:
    :param out_dir:
    :return:
    """

    # Paramters
    neg_ctrls = ("DMSO",)
    neg_ctrls_concentrations = (None,)
    pos_ctrls = ("Bortezomib", "Irinotecan / SN-38",
                 "Volasertib", "Methotrexate", "Staurosporine_500nM")
    pos_ctrls_concentrations = ((0.2, 1), (0.2, 1), (0.2, 1),
                                (0.2, 1), None)

    # Check if training data already exists
    hdf5_fn = os.path.join(out_dir, "TrainingData.h5")
    if os.path.isfile(hdf5_fn):
        return False

    proj_dir = os.path.join(promise_dir, "hdf5projection")
    all_cell_lines = sorted(set([
        s[0:7] for s in os.listdir(proj_dir) if
        s.startswith("D0")]))

    # Load wells
    lib_fn = os.path.join(lib_dir, "python_friendly", "L08.xlsx")
    layout = pd.read_excel(lib_fn)
    neg_ctrls_wells = []
    for neg_ctrl_i in range(len(neg_ctrls)):
        neg_ctrl = neg_ctrls[neg_ctrl_i]
        neg_ctrl_concentrations = neg_ctrls_concentrations[neg_ctrl_i]
        if neg_ctrl_concentrations is None:
            neg_ctrls_wells.append(
                layout.loc[layout["Product.Name"].values.astype(str) == neg_ctrl,
                           "Well_ID_384"].values.astype(str))
        else:
            neg_ctrls_wells.append(layout.loc[
                np.logical_and(
                    layout["Product.Name"].values.astype(str) == neg_ctrl,
                    layout["concentration"].values.astype(str) == neg_ctrl_concentrations),
                "Well_ID_384"].values.astype(str))
    neg_ctrls_wells = np.concatenate(neg_ctrls_wells)
    pos_ctrls_wells = []
    for pos_ctrl_i in range(len(pos_ctrls)):
        pos_ctrl = pos_ctrls[pos_ctrl_i]
        pos_ctrl_concentrations = pos_ctrls_concentrations[pos_ctrl_i]
        if pos_ctrl_concentrations is None:
            pos_ctrls_wells.append(
                layout.loc[layout["Product.Name"].values.astype(str) == pos_ctrl,
                           "Well_ID_384"].values.astype(str))
        else:
            pos_ctrls_wells.append(layout.loc[
                np.logical_and(
                    layout["Product.Name"].values.astype(str) == pos_ctrl,
                    np.in1d(layout["concentration"].values.astype(float),
                            np.array(pos_ctrl_concentrations))),
                "Well_ID_384"].values.astype(str))
    pos_ctrls_wells = np.concatenate(pos_ctrls_wells)

    # Match the number of wells to use
    max_len = min(len(neg_ctrls_wells), len(pos_ctrls_wells))
    pos_ctrls_wells = random.sample(pos_ctrls_wells, max_len)
    neg_ctrls_wells = random.sample(neg_ctrls_wells, max_len)

    pos_patches = []
    neg_patches = []
    for cell_line in all_cell_lines:
        l08_plates = [
            s for s in os.listdir(proj_dir)
            if s.startswith(cell_line)
            if s.endswith("L08")
            if os.path.isdir(os.path.join(proj_dir, s))]

        for plate in l08_plates:
            for pos_well in pos_ctrls_wells:
                print plate, pos_well
                pos_patches.append(create_patches(
                    plate, pos_well[0], pos_well[1:3],
                    patch_rad, patches_per_well, promise_dir))
            for neg_well in neg_ctrls_wells:
                print plate, neg_well
                neg_patches.append(create_patches(
                    plate, neg_well[0], neg_well[1:3],
                    patch_rad, patches_per_well, promise_dir))

    pos_patches = np.concatenate(pos_patches)
    neg_patches = np.concatenate(neg_patches)

    # Match sample sizes
    max_len = min(len(pos_patches), len(neg_patches))
    pos_patches = pos_patches[np.random.choice(
        range(len(pos_patches)), max_len, replace=False), ...]
    neg_patches = neg_patches[np.random.choice(
        range(len(neg_patches)), max_len, replace=False), ...]

    # Set up training data
    x = np.concatenate((pos_patches, neg_patches), axis=0)
    y = np.concatenate((
        np.array([[1, 0]] * len(pos_patches)),
        np.array([[0, 1]] * len(neg_patches))))

    with h5py.File(hdf5_fn, "w-") as h5handle:
        h5handle.create_dataset(name="images", data=x)
        h5handle.create_dataset(name="labels", data=y)

    return True
