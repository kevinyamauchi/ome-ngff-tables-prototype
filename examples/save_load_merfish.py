import napari
import re
import anndata as ad
import json
import os
import shutil
from ngff_tables_prototype.writer import write_spatial_anndata
from ngff_tables_prototype.reader import load_to_napari_viewer
import numpy as np

output_fpath = "test_merfish.zarr"


def write_merfish_adata() -> None:
    data_folder = "../notebooks/merfish/output"
    files = {
        "cells": "cells.h5ad",
        "image": "image.npy",
        "image_transform": "image_transform.json",
        "points": "points.h5ad",
        "polygons": "polygons.h5ad",
    }

    def pj(file):
        return os.path.join(data_folder, file)

    if not all([os.path.isfile(pj(file)) for file in files.values()]):
        raise FileNotFoundError(
            "run merfish.ipynb or merfish.py to generate the required files"
        )

    if True:
        if os.path.isdir(output_fpath):
            shutil.rmtree(output_fpath)

    # load the data
    # image
    image = np.load(pj(files["image"]))
    j = json.load(open(pj(files["image_transform"]), "r"))
    image_translation = np.array([j["translation_x"], j["translation_y"]])
    image_scale_factors = np.array([j["scale_factor_x"], j["scale_factor_y"]])

    # cells
    a_cells = ad.read_h5ad(pj(files["cells"]))

    # points
    a_points = ad.read_h5ad(pj(files["points"]))

    # polygons
    a_polygons = ad.read(pj(files["polygons"]))

    if not os.path.isdir(output_fpath):
        write_spatial_anndata(
            file_path=output_fpath,
            # image
            image_axes=["y", "x"],
            image=image,
            image_translation=image_translation,
            image_scale_factors=image_scale_factors,
            # expression
            tables_adata=a_cells,
            tables_region="circles/circles_table",
            # circles
            circles_adata=a_cells,
            # points
            points_adata=a_points,
            # polygons
            polygons_adata=a_polygons,
        )

    return


if __name__ == "__main__":
    write_merfish_adata()
    viewer = load_to_napari_viewer(
        file_path=output_fpath,
        # groups=["circles/circles_table", 'polygons/polygons_table'],
        groups=[
            "circles/circles_table",
            "tables/regions_table",
            "points/points_table",
            "polygons/polygons_table",
        ],
    )
    napari.run()
