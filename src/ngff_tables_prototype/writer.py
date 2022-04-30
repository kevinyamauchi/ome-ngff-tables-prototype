from typing import Optional, Union, Tuple, Any, List, Dict

import anndata
from anndata import AnnData
from anndata._io.utils import write_attribute
import numpy as np
import pandas as pd
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image, write_labels
from scipy import sparse
import zarr

from .table_utils import df_to_anndata


def write_table(
        group: zarr.Group,
        adata: AnnData,
        chunks=None,
        **dataset_kwargs,
):
    """code from: https://github.com/theislab/anndata/blob/0.7.8/anndata/_io/zarr.py
    """
    # AnnData requires an upper case key for X
    # todo: make context manager?
    group.store.normalize_keys = False

    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    # TODO: Use spec writing system for this
    group.attrs.setdefault("encoding-type", "anndata")
    group.attrs.setdefault("encoding-version", "0.1.0")
    if chunks is not None and not isinstance(adata.X, sparse.spmatrix):
        write_attribute(group, "X", adata.X, dict(chunks=chunks, **dataset_kwargs))
    else:
        write_attribute(group, "X", adata.X, dataset_kwargs)
    write_attribute(group, "obs", adata.obs, dataset_kwargs)
    write_attribute(group, "var", adata.var, dataset_kwargs)
    write_attribute(group, "obsm", adata.obsm, dataset_kwargs)
    write_attribute(group, "varm", adata.varm, dataset_kwargs)
    write_attribute(group, "obsp", adata.obsp, dataset_kwargs)
    write_attribute(group, "varp", adata.varp, dataset_kwargs)
    write_attribute(group, "layers", adata.layers, dataset_kwargs)
    write_attribute(group, "uns", adata.uns, dataset_kwargs)
    write_attribute(group, "raw", adata.raw, dataset_kwargs)

    # revert the store to normalizing keys
    group.store.normalize_keys = True


def write_points_dataset(
    file_path: str,
    image: Optional[np.ndarray] = None,
    image_chunks: Union[Tuple[Any, ...], int] = None,
    image_axes: Union[str, List[str]] = None,
    points: Optional[pd.DataFrame] = None,
    points_dense_columns: Optional[List[str]] = None,
    points_var: Optional[pd.DataFrame] = None,
    points_group_name: str = 'points'
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
        write_image(
            image=image,
            group=root,
            chunks=image_chunks,
            axes=image_axes
        )

    if points is not None:
        points_group = root.create_group(points_group_name)
        points_anndata = df_to_anndata(
            df=points,
            dense_columns=points_dense_columns,
            var=points_var
        )
        write_table(points_group, points_anndata)

        # add a flag to the group attrs to denote this is points data
        points_group.attrs['@type'] = 'ngff:points'


def write_spatial_anndata(
    file_path: str,
    image: Optional[np.ndarray]=None,
    image_chunks: Union[Tuple[Any, ...], int] = None,
    image_axes: Union[str, List[str]] = None,
    label_image: Optional[np.ndarray]=None,
    label_group_name: str = "labels",
    label_name: str = "labels",
    label_column : Optional[str] = "cell_id",
    adata: Optional[anndata.AnnData] =None,
    table_group_name: str = "table",
):
    """ Write a spatial anndata object to ome-zarr
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
    """
    # create the zarr root
    store = parse_url(file_path, mode="w").store
    root = zarr.group(store=store)

    if image is not None:
        # if necessary add the image
        write_image(
            image=image,
            group=root,
            chunks=image_chunks,
            axes=image_axes
        )

    if label_image is not None:
        # label_image_group = root.create_group(label_group_name)
        write_labels(
            label_image,
            group=root,
            name=label_name
        )

    if adata is not None:
        from anndata.experimental import write_elem

        # write_table(table_group, adata)
        write_elem(root, table_group_name, adata)
        table_group = root[table_group_name]
        table_group.attrs["@type"] = "ngff:table"
        table_group.attrs["label"] = label_group_name
        table_group.attrs["label_column"] = label_column
