import os

from anndata import AnnData
from anndata._io import read_zarr
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
import napari
from napari_ome_zarr._reader import transform
from napari.types import LayerDataTuple


def load_points_to_anndata(file_path: str, points_group: str) -> AnnData:
    return read_zarr(os.path.join(file_path, points_group))


def _anndata_to_napari_points(anndata_obj: AnnData) -> LayerDataTuple:
    points_coords = anndata_obj.X
    layer_properties = anndata_obj.obs
    return points_coords, {'properties': layer_properties}, 'points'


def load_to_napari_layer_data(file_path: str, points_group: str):
    ome_zarr = parse_url(file_path)
    reader = Reader(ome_zarr)
    reader_func = transform(reader())
    layer_data = reader_func(file_path)

    # load points points
    anndata_obj = load_points_to_anndata(file_path, points_group)
    layer_data.append(_anndata_to_napari_points(anndata_obj))

    return layer_data

def load_to_napari_viewer(file_path: str, points_group: str) -> napari.Viewer:
    layer_data = load_to_napari_layer_data(file_path, points_group)

    viewer = napari.Viewer()

    for layer in layer_data:
        viewer._add_layer_from_data(*layer)
    return viewer
