import numpy as np
import pandas as pd

from ngff_tables_prototype.writer import write_points_dataset

# create the image
rng = np.random.default_rng(42)
image = rng.poisson(64, size=(128, 128)).astype(np.uint8)

# create points data
n_points = 10
n_dim = 2
points_coords = 100 * rng.random((n_points, n_dim), dtype=float)
points = pd.DataFrame(
    {
        'y': points_coords[:, 0],
        'x': points_coords[:, 1],
        'size': np.repeat(3, n_points)
    }
)
var = pd.DataFrame({'axis': ['y', 'x'], 'scale': [1, 1]})

print(points)

output_fpath = 'test_image'
write_points_dataset(
    file_path=output_fpath,
    image=image,
    image_axes='yx',
    points=points,
    points_dense_columns=['y', 'x'],
    points_var=var
)
