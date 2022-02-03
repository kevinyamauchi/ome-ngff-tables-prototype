import os

import napari
import numpy as np
import pandas as pd

from ngff_tables_prototype.writer import write_points_dataset
from ngff_tables_prototype.reader import load_points_to_anndata, load_to_napari_viewer

# create the image
rng = np.random.default_rng(42)
image = rng.poisson(64, size=(128, 128)).astype(np.uint8)

# create points data
n_points = 10
n_dim = 2
points_coords = 100 * rng.random((n_points, n_dim), dtype=np.float32)
points = pd.DataFrame(
    {
        'y': points_coords[:, 0],
        'x': points_coords[:, 1],
        'size': np.repeat(3, n_points)
    }
)
var = pd.DataFrame({'axis': ['y', 'x'], 'scale': [1, 1]})

print(points)

output_fpath = 'test_image.zarr'

if not os.path.isdir(output_fpath):
    write_points_dataset(
        file_path=output_fpath,
        image=image,
        image_axes='yx',
        points=points,
        points_dense_columns=['y', 'x'],
        points_var=var
    )

# reload the points object to anndata and confirm the data are the same
anndata_obj = load_points_to_anndata(file_path='test_image.zarr', points_group='points')
np.testing.assert_almost_equal(points_coords, anndata_obj.X)

viewer = load_to_napari_viewer(file_path='test_image.zarr', points_group='points')
napari.run()
