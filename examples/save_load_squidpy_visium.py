import os
import shutil
import napari
import squidpy as sq
from ngff_tables_prototype.writer import write_spatial_anndata
from ngff_tables_prototype.reader import load_to_napari_viewer
import numpy as np
from anndata import AnnData

output_fpath = "test_visium.zarr"


def write_visium_adata() -> None:
    adata = sq.datasets.visium_fluo_adata_crop()

    lib_id = "V1_Adult_Mouse_Brain_Coronal_Section_2"
    spatial_key = "spatial"

    image = adata.uns[spatial_key][lib_id]["images"]["hires"]
    scalef = adata.uns[spatial_key][lib_id]["scalefactors"]["tissue_hires_scalef"]
    adata.obsm[spatial_key] = adata.obsm[spatial_key] * scalef

    adata.X = adata.X.A.copy()

    tables_adata = adata.copy()
    circles_adata = AnnData(
        None,
        obs=adata.obs.copy(),
        var=adata.var.copy(),
        obsm={"spatial": adata.obsm["spatial"]},
        uns=adata.uns.copy(),
    )
    circles_adata.obs_names = adata.obs_names.copy()
    circles_adata.var_names = adata.var_names.copy()

    if True:
        if os.path.isdir(output_fpath):
            shutil.rmtree(output_fpath)

    if not os.path.isdir(output_fpath):
        write_spatial_anndata(
            file_path=output_fpath,
            image_axes=["c", "y", "x"],
            image=np.swapaxes(image, 2, 0),
            circles_adata=circles_adata,
            tables_adata=tables_adata,
            tables_region="circles/circles_table",
        )

    return


if __name__ == "__main__":
    write_visium_adata()
    viewer = load_to_napari_viewer(
        file_path=output_fpath,
        groups=["circles/circles_table", "tables/regions_table"],
    )
    napari.run()
