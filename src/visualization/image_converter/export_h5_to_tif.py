from h5py import File
from imageio import imwrite
import os
import sys
import subprocess

h5_remote_path = sys.argv[1]

h5_file_name = os.path.basename(h5_remote_path)

subprocess.run(["scp", "remote_1@b110-imac31:/Volumes/b110_data/B110-Isilon2"+h5_remote_path, h5_file_name])

prefix = h5_file_name.split(".")[0]+"_TIFF/"

if not os.path.exists(prefix):
    os.makedirs(prefix)

h5_file = File(h5_file_name, 'r')

image_set = h5_file['images']

for field in range(4):
    for channel in range(2):
        imwrite('{}{}_export_fld{}_ch{}.tif'.format(prefix, h5_file_name.split(".")[0], field + 1, channel + 1),
                image_set[field, channel])
