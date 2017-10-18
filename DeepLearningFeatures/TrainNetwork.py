# This script trains a network to tell the difference between positive and
# negative patches
from __future__ import division
import random
import pickle
import os
import numpy as np
import pandas as pd
import h5py
import ImageMosaic
import sys
import skimage.morphology
import skimage.filters.rank
import sklearn.model_selection
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
import keras.backend
sys.path.append("/Users/jansauer/Thesis/Projects/PROMISE/DeepLearningFeatures")
import ImagePreprocessor


def create_patches(image_fn, patch_rad=40, num_patches=10000):
    """
    Load an image and select patches centered on high-contrast areas
    :param image_fn:
    :param patch_rad:
    :param num_patches:
    :return:
    """
    with h5py.File(image_fn, "r") as h5h:
        _proj = h5h["images"][()]
    _proj = _proj.transpose((0, 2, 3, 1))
    _proj = ImageMosaic.create_mosaic(_proj, border_size=0)
    _proj = ImagePreprocessor.preprocess_image(image=_proj)
    _proj_gs = np.mean(_proj, axis=-1)
    _proj_q = np.percentile(_proj_gs, 99.9)
    _proj_gs[_proj_gs > _proj_q] = _proj_q
    _proj_gs -= _proj_gs.min()
    _proj_gs /= _proj_gs.max()
    _proj_entropy = skimage.filters.rank.entropy(
        image=_proj_gs, selem=skimage.morphology.disk(15))
    _proj_entropy_mask = _proj_entropy >= np.percentile(_proj_entropy, 80)
    _points = np.transpose(np.where(_proj_entropy_mask))
    _points = np.stack(random.sample(_points, num_patches))
    _points = _points[np.logical_and(
        patch_rad <= _points[..., 0],
        _points[..., 0] < (_proj.shape[0] - patch_rad)), ...]
    _points = _points[np.logical_and(
        patch_rad <= _points[..., 1],
        _points[..., 1] < (_proj.shape[1] - patch_rad)), ...]
    _patches = []
    for pp in _points:
        _patches.append(_proj[
                           (pp[0] - patch_rad):(pp[0] + patch_rad + 1),
                           (pp[1] - patch_rad):(pp[1] + patch_rad + 1)])
    return np.stack(_patches)


def get_wells_from_layout(layout_file, treatments, concentrations=None):
    """
    Load control plates from library file
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


def create_training_data(
        target_patch_numbers=50000,
        target_patch_radius=80,
        proj_dir="/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5projection",
        lib_dir="/collab-ag-fischer/PROMISE/layouts", neg_ctrls=("DMSO",),
        pos_ctrls=("Bortezomib", "Irinotecan / SN-38",
                   "Volasertib", "Methotrexate", "Staurosporine_500nM"),
        out_dir="."):
    """
    Create a training data hdf5 file from the libraries. This currently only
    uses the L08 plates, as only this library has positive and negative
    controls in known concentrations.
    :param proj_dir:
    :param lib_dir:
    :param neg_ctrls:
    :param pos_ctrls:
    :return:
    """

    hdf5_fn = os.path.join(out_dir, "TrainingData.h5")
    if os.path.isfile(hdf5_fn):
        return False

    all_cell_lines = sorted(set([
        s[0:7] for s in os.listdir(proj_dir) if
        s.startswith("D0")]))

    lib_fn = os.path.join(lib_dir, "python_friendly", "L08.xlsx")
    layout = pd.read_excel(lib_fn)
    neg_ctrl = layout.loc[
        np.in1d(layout["Product.Name"].values.astype(str),
                np.array(neg_ctrls)), "Well_ID_384"].values.astype(str)
    pos_ctrl = layout.loc[np.logical_and(
        np.in1d(layout["Product.Name"].values.astype(str),
                np.array(pos_ctrls)),
        np.in1d(layout["concentration"].values.astype(float),
                (0.2, 1))), "Well_ID_384"].values.astype(str)

    # Match the number of wells to use
    max_len = min(len(neg_ctrl), len(pos_ctrl))
    pos_ctrl = random.sample(pos_ctrl, max_len)
    neg_ctrl = random.sample(neg_ctrl, max_len)

    pos_patches = []
    neg_patches = []
    for cell_line in all_cell_lines[0:1]:
        l08_plates = [
            s for s in os.listdir(proj_dir)
            if s.startswith(cell_line)
            if s.endswith("L08")
            if os.path.isdir(os.path.join(proj_dir, s))]

        num_patches_per_well = int(np.ceil(
            target_patch_numbers /
            (max_len * len(all_cell_lines) * len(l08_plates))))

        for plate in l08_plates:
            for pos_well in pos_ctrl:
                print cell_line, plate, pos_well
                pos_well_fn = os.path.join(
                    proj_dir, plate,
                    "%s_%s_%s_contrastProjections.h5" %
                    (plate, pos_well[0], pos_well[1:3]))
                pos_patches.append(create_patches(
                    image_fn=pos_well_fn, num_patches=num_patches_per_well))
            for neg_well in neg_ctrl:
                print cell_line, plate, neg_well
                neg_well_fn = os.path.join(
                    proj_dir, plate,
                    "%s_%s_%s_contrastProjections.h5" %
                    (plate, neg_well[0], neg_well[1:3]))
                neg_patches.append(create_patches(
                    image_fn=neg_well_fn, num_patches=num_patches_per_well))

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


def train_network(out_dir="."):
    hdf5_fn = os.path.join(out_dir, "TrainingData.h5")
    if not os.path.isfile(hdf5_fn):
        create_training_data()
    with h5py.File(hdf5_fn, "r") as h5handle:
        x = h5handle["images"][()]
        y = h5handle["labels"][()]

    # Split into training and testing sets
    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(
        x, y, test_size=0.25)

    # Set up and run network
    model = Sequential()
    model.add(Conv2D(32, kernel_size=(3, 3),
                     activation='relu',
                     input_shape=x_train.shape[1:4]))
    model.add(Conv2D(64, (3, 3), activation='relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))
    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(y.shape[1], activation='softmax'))

    model.compile(loss=keras.losses.categorical_crossentropy,
                  optimizer=keras.optimizers.Adadelta(),
                  metrics=['accuracy'])

    model.fit(x_train, y_train,
              batch_size=128,
              epochs=1,
              verbose=1,
              validation_data=(x_test, y_test))
    model.save(os.path.join(out_dir, "KerasModel.h5"))
    score = model.evaluate(x_test, y_test, verbose=0)
    with open(os.path.join(out_dir, "KerasModelAcc.pkl"), "w-") as f:
        pickle.dump(score, f, pickle.HIGHEST_PROTOCOL)
    # print('Test loss:', score[0])
    # print('Test accuracy:', score[1])


def apply_network(out_dir="."):
    """
    Apply the trained network to the positive and negative controls
    :return:
    """
    proj_dir = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5projection"
    lib_dir = "/collab-ag-fischer/PROMISE/layouts"

    # Load network
    model = keras.models.load_model(os.path.join(out_dir, "KerasModel.h5"))

    # Prep extraction function
    get_features = keras.backend.function(
        [model.layers[0].input, keras.backend.learning_phase()],
        [model.layers[-2].output])

    # Get wells for controls
    cell_line = "D004T01"
    l08_plates = [
        s for s in os.listdir(proj_dir)
        if s.startswith(cell_line)
        if s.endswith("L08")
        if os.path.isdir(os.path.join(proj_dir, s))]

    layout_fn = os.path.join(lib_dir, "python_friendly", "L08.xlsx")
    neg_ctrls = ("DMSO",)
    neg_ctrl_wells = [
        p + "_" + w for p in l08_plates for w in
        get_wells_from_layout(layout_fn, neg_ctrls)]
    pos_ctrls = ("Bortezomib", "Irinotecan / SN-38",
                 "Volasertib", "Methotrexate", "Staurosporine_500nM")
    pos_ctrl_wells = [
        p + "_" + w for p in l08_plates for w in
        get_wells_from_layout(layout_fn, pos_ctrls)]

    features_pos = []
    for pcw in pos_ctrl_wells[0:2]:
        plate, well = pcw.split("_")
        image_fn = os.path.join(
            proj_dir, plate,
            "%s_%s_%s_contrastProjections.h5" %
            (plate, well[0], well[1:3]))

        patches = create_patches(image_fn, num_patches=100)
        features_pos.append(get_features([patches, 0])[0])
    features_pos = np.concatenate(features_pos, axis=0)
    np.save(os.path.join(out_dir, "Features_PosCtrl.npy"), features_pos)

    features_neg = []
    for ncw in neg_ctrl_wells:
        plate, well = ncw.split("_")
        image_fn = os.path.join(
            proj_dir, plate,
            "%s_%s_%s_contrastProjections.h5" %
            (plate, well[0], well[1:3]))

        patches = create_patches(image_fn, num_patches=100)
        features_neg.append(get_features([patches, 0]))
    features_neg = np.concatenate(features_neg, axis=0)
    np.save(os.path.join(out_dir, "Features_NegCtrl.npy"), features_neg)