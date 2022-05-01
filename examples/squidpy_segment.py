import squidpy as sq
from ngff_tables_prototype.writer import write_spatial_anndata
import numpy as np


def write_segment_adata() -> None:
    adata = sq.datasets.mibitof()
    lib_id = "point8"
    spatial_key = "spatial"
    adata = adata[adata.obs.library_id == lib_id].copy()
    image = adata.uns[spatial_key][lib_id]["images"]["hires"]
    segment = adata.uns[spatial_key][lib_id]["images"]["segmentation"]
    adata.X = adata.X.A.copy()
    tables_adata = adata.copy()

    write_spatial_anndata(
        file_path="test_segment.zarr",
        image_axes=["c", "y", "x"],
        image=np.swapaxes(image, 2, 0),
        label_image=segment,
        tables_adata=tables_adata,
        tables_instance_key="cell_id",
    )

    return


if __name__ == "__main__":
    write_segment_adata()
