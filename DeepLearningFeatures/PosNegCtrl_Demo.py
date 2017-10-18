# This is a demo script to train a neural network to differentiate between
# positive and negative controls. To permit generality, I want to use no
# segmentation results for this.
#
# It turns out that Bortezomib at the highest concentration is a far
# better positive control than Staurosporine based on organoid count
# and size.
#
# This demo script uses only human lines and library L08.

from __future__ import division
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import skimage.filters.rank
import skimage.morphology
import ImageMosaic
import random
import h5py
import sklearn.model_selection
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D
import keras.backend
sys.path.append("/Users/jansauer/Thesis/Projects/PROMISE/DeepLearningFeatures")
import ImagePreprocessor


# Environment variables
PROJDIR = "/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5projection"
LIBDIR = "/collab-ag-fischer/PROMISE/layouts"
OUTDIR = "/Users/jansauer/Thesis/Projects/PROMISE/DeepLearningFeatures"


def create_patches(image_fn, patch_rad=40, num_patches=10000):
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
        _points[..., 0] < (_proj.shape[0] - patch_rad)),]
    _points = _points[np.logical_and(
        patch_rad <= _points[..., 1],
        _points[..., 1] < (_proj.shape[1] - patch_rad)),]
    _patches = []
    for pp in _points:
        _patches.append(_proj[
                           (pp[0] - patch_rad):(pp[0] + patch_rad + 1),
                           (pp[1] - patch_rad):(pp[1] + patch_rad + 1)])
    return np.stack(_patches)


def create_training_data():
    all_cell_lines = sorted(set([
        s[0:7] for s in os.listdir(PROJDIR) if
        s.startswith("D0")]))

    # Only L08 has both positive and negative controls so only these plates are used
    cell_line = all_cell_lines[0]
    l08_plates = [
        s for s in os.listdir(PROJDIR)
        if s.startswith(cell_line)
        if s.endswith("L08")
        if os.path.isdir(os.path.join(PROJDIR, s))]

    lib_fn = os.path.join(LIBDIR, "python_friendly", "L08.xlsx")
    layout = pd.read_excel(lib_fn)
    neg_ctrl = layout.loc[
        layout["Product.Name"] == "DMSO", "Well_ID_384"].values.astype(str)
    pos_ctrl = layout.loc[
        (layout["Product.Name"] == "Bortezomib") & (layout["concentration"] == 1),
        "Well_ID_384"].values.astype(str)

    # Match the number of wells to use
    max_len = min(len(neg_ctrl), len(pos_ctrl))
    pos_ctrl = random.sample(pos_ctrl, max_len)
    neg_ctrl = random.sample(neg_ctrl, max_len)

    # Get positive and negative patches
    pos_patches = []
    neg_patches = []
    for plate in l08_plates:
        for pos_well in pos_ctrl:
            pos_well_fn = os.path.join(
                PROJDIR, plate,
                "%s_%s_%s_contrastProjections.h5" %
                (plate, pos_well[0], pos_well[1:3]))
            pos_patches.append(create_patches(image_fn=pos_well_fn))
        for neg_well in neg_ctrl:
            neg_well_fn = os.path.join(
                PROJDIR, plate,
                "%s_%s_%s_contrastProjections.h5" %
                (plate, neg_well[0], neg_well[1:3]))
            neg_patches.append(create_patches(image_fn=neg_well_fn))

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
    return x, y


hdf5_fn = os.path.join(OUTDIR, "TrainingData.h5")
if os.path.isfile(hdf5_fn):
    with h5py.File(hdf5_fn, "r") as h5handle:
        x = h5handle["images"][()]
        y = h5handle["labels"][()]
else:
    x, y = create_training_data()
    with h5py.File(hdf5_fn, "w-") as h5handle:
        h5handle.create_dataset(name="images", data=x)
        h5handle.create_dataset(name="labels", data=y)

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
score = model.evaluate(x_test, y_test, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])


# Create a function to get FC outputs
get_fc_outputs = keras.backend.function(
    inputs=[model.layers[0].input, keras.backend.learning_phase()],
    outputs=[model.layers[-3].output])

fc_output = get_fc_outputs([x[y[..., 0] == 1, ...]])


