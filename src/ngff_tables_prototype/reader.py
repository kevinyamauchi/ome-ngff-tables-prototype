import os

from anndata import AnnData
from anndata._io import read_zarr
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
import napari
from napari_ome_zarr._reader import transform
from napari.types import LayerDataTuple
import pandas as pd
import zarr


def load_table_to_anndata(file_path: str, table_group: str) -> AnnData:
    return read_zarr(os.path.join(file_path, table_group))


def _anndata_to_napari_points(anndata_obj: AnnData) -> LayerDataTuple:
    points_coords = anndata_obj.X
    layer_properties = anndata_obj.obs
    return points_coords, {'properties': layer_properties}, 'points'

def _anndata_to_napari_labels_features(anndata_obj: AnnData) -> pd.DataFrame:
    df = pd.concat([anndata_obj.to_df(), anndata_obj.obs], axis=1)

    # hacky add dummy data to first row to meet the condition for napari features
    # https://github.com/napari/napari/blob/2648c307168fdfea7b5c03aaad4590972e5dfe2c/napari/layers/labels/labels.py#L72-L74
    first_row = df.iloc[0, :]
    first_row["cell_id"] = 0

    return pd.concat([first_row, df])

def load_to_napari_layer_data(file_path: str, table_group_key: str):
    ome_zarr = parse_url(file_path)
    reader = Reader(ome_zarr)
    reader_func = transform(reader())
    layer_data = reader_func(file_path)

    # load table
    anndata_obj = load_table_to_anndata(file_path, table_group_key)
    table_group = zarr.group(ome_zarr.store)[table_group_key]
    if "label" in table_group.attrs:
        features_table = _anndata_to_napari_labels_features(anndata_obj)

        # get the label item
        label_name = table_group.attrs["label"]
        for i, (data, kwargs, layer_type) in enumerate(layer_data):

            if (layer_type.lower() == "labels") and (kwargs["name"] == label_name):
                break
        kwargs.update({"features": features_table})
        layer_data[i] = (data, kwargs, layer_type)

    return layer_data

def load_to_napari_viewer(file_path: str, points_group: str) -> napari.Viewer:
    layer_data = load_to_napari_layer_data(file_path, points_group)

    viewer = napari.Viewer()

    for layer in layer_data:
        viewer._add_layer_from_data(*layer)
    return viewer
