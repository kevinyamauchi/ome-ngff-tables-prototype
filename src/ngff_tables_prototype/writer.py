from typing import Optional, Union, Tuple, Any, List

import anndata
from anndata import AnnData
import numpy as np
import pandas as pd
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image, write_labels
import zarr
from anndata.experimental import write_elem
from .table_utils import df_to_anndata

# make function that writes region table
# function that writes region circles


def write_table_points(
    group: zarr.Group,
    adata: AnnData,
    table_group_name: str = "points_table",
    group_type: str = "ngff:points_table",
):
    write_elem(group, table_group_name, adata)
    table_group = group[table_group_name]
    table_group.attrs["@type"] = group_type


def write_table_polygons(
    group: zarr.Group,
    adata: AnnData,
    table_group_name: str = "shapes_table",
    group_type: str = "ngff:shapes_table",
):
    pass


def write_table_circles(
    group: zarr.Group,
    adata: AnnData,
    table_group_name: str = "circles_table",
    group_type: str = "ngff:circles_table",
):
    write_elem(group, table_group_name, adata)
    table_group = group[table_group_name]
    table_group.attrs["@type"] = group_type


def write_table_regions(
    group: zarr.Group,
    adata: AnnData,
    table_group_name: str = "regions_table",
    group_type: str = "ngff:regions_table",
    region: Union[str, List[str]] = "features",
    region_key: Optional[str] = None,
    instance_key: Optional[str] = None,
):
    write_elem(group, table_group_name, adata)
    table_group = group[table_group_name]
    table_group.attrs["@type"] = group_type
    table_group.attrs["region"] = region
    table_group.attrs["region_key"] = region_key
    table_group.attrs["instance_key"] = instance_key


def write_points_dataset(
    file_path: str,
    image: Optional[np.ndarray] = None,
    image_chunks: Union[Tuple[Any, ...], int] = None,
    image_axes: Union[str, List[str]] = None,
    points: Optional[pd.DataFrame] = None,
    points_dense_columns: Optional[List[str]] = None,
    points_var: Optional[pd.DataFrame] = None,
    points_group_name: str = "points",
):
    """Write a dataset with image and points to ome-zarr

    Parameters
    ----------
    file_path : str
        path to save the zarr file
    image : Optional[np.ndarray]
        image array to save. if None, no image is saved.
        Default value is None.
    image_chunks : Union[Tuple[Any, ...], int]
        Chunking for the image data. See ome-zarr-py for details.
    image_axes : Union[str, List[str]]
        The labels for the image axes. See ome-zarr-py for details.
    points : Optional[pd.DataFrame]
        Table of points to store.
    points_dense_columns : List[str]
        THe names of the columns in points to store in the dense array.
        Must be the names of the coordinate columns.
    points_var : Optional[pd.DataFrame]
        The annotations for the columns of the dense matrix X that
        get stored in AnnData.var. Must be ordered to match
        dense_columns.
    points_group_name : str
        The name of the group to store the points in.
        Default value is 'points'
    """
    # create the zarr root
    store = parse_url(file_path, mode="w").store
    root = zarr.group(store=store)

    if image is not None:
        # if necessary add the image
        write_image(image=image, group=root, chunks=image_chunks, axes=image_axes)

    if points is not None:
        points_group = root.create_group(points_group_name)
        points_anndata = df_to_anndata(
            df=points, dense_columns=points_dense_columns, var=points_var
        )
        write_table(points_group, points_anndata)

        # add a flag to the group attrs to denote this is points data
        points_group.attrs["@type"] = "ngff:points"


def write_spatial_anndata(
    file_path: str,
    # image group
    image: Optional[np.ndarray] = None,
    image_chunks: Union[Tuple[Any, ...], int] = None,
    image_axes: Union[str, List[str]] = None,
    # label group
    label_image: Optional[np.ndarray] = None,
    # label_name: str = "label_image",
    # table group
    tables_adata: Optional[AnnData] = None,
    tables_region: Optional[Union[str, List[str]]] = None,
    tables_region_key: Optional[str] = None,
    tables_instance_key: Optional[str] = None,
    # shape group
    circles_adata: Optional[AnnData] = None,
    polygons_adata: Optional[AnnData] = None,
    # points group
    points_adata: Optional[AnnData] = None,
):
    """Write a spatial anndata object to ome-zarr
    Parameters
    ----------
    file_path : str
        path to save the zarr file
    image : Optional[np.ndarray]
        image array to save. if None, no image is saved.
        Default value is None.
    image_chunks : Union[Tuple[Any, ...], int]
        Chunking for the image data. See ome-zarr-py for details.
    image_axes : Union[str, List[str]]
        The labels for the image axes. See ome-zarr-py for details.
    label_image : Union[str, List[str]]
        The label image (raster-mask). See ome-zarr-py for details.
    tables_adata:
        The :class:`anndata.AnnData` table with gene expression and annotations.
    tables_region
        The :class:`anndata.AnnData` region table that maps to one (or more) label image.
    tables_region_key
        The key in :attr:`AnnData.obs` that stores unique pointers to one (or more) label image.
    tables_instance_key
        The key in :attr:`AnnData.obs` that stores unique pointers to regions.
    circles_adata
        The :class:`anndata.AnnData` circle table that store coordinates of circles and metadata.
    polygons_adata
        The :class:`anndata.AnnData` polygons table that store coordinates of polygons and metadata.
    points_adata
        The :class:`anndata.AnnData` point table that store coordinates of points and metadata.

    Returns
    -------
    Nothing, save the file in Zarr store.
    """
    # create the zarr root
    store = parse_url(file_path, mode="w").store
    root = zarr.group(store=store)

    if image is not None:
        # if necessary add the image
        write_image(image=image, group=root, chunks=image_chunks, axes=image_axes)

    if label_image is not None:
        # label_image_group = root.create_group(label_group_name)
        write_labels(label_image, group=root, name=label_name)

    if tables_adata is not None:
        write_table_regions(
            tables_adata,
            group=root,
            tables_region=tables_region,
            tables_region_key=tables_region_key,
            tables_instance_key=tables_instance_key,
        )
    if circles_adata is not None:
        write_table_regions(
            circles_adata,
            group=root,
        )
    if polygons_adata is not None:
        write_table_regions(
            polygons_adata,
            group=root,
        )
    if points_adata is not None:
        write_table_regions(
            points_adata,
            group=root,
        )
