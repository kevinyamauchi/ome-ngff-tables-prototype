import napari
import os
import shutil
import squidpy as sq
from ngff_tables_prototype.writer import write_spatial_anndata
from ngff_tables_prototype.reader import load_to_napari_viewer
import numpy as np

output_fpath = "test_segment.zarr"


def write_segmentation_adata() -> None:
    adata = sq.datasets.mibitof()
    lib_id = "point8"
    spatial_key = "spatial"
    adata = adata[adata.obs.library_id == lib_id].copy()
    image = adata.uns[spatial_key][lib_id]["images"]["hires"]
    segment = adata.uns[spatial_key][lib_id]["images"]["segmentation"]
    adata.X = adata.X.A.copy()
    tables_adata = adata.copy()

    if True:
        if os.path.isdir(output_fpath):
            shutil.rmtree(output_fpath)

    if not os.path.isdir(output_fpath):
        write_spatial_anndata(
            file_path=output_fpath,
            image_axes=["c", "y", "x"],
            image=np.swapaxes(image, 2, 0),
            label_image=segment,
            tables_adata=tables_adata,
            tables_region="labels/label_image",
            tables_instance_key="cell_id",
        )

    return


if __name__ == "__main__":
    write_segmentation_adata()
    viewer = load_to_napari_viewer(
        file_path=output_fpath,
        groups=["labels/label_image", "tables/regions_table"],
    )
    napari.run()
