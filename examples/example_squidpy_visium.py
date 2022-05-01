import squidpy as sq
from ngff_tables_prototype.writer import write_spatial_anndata
import numpy as np
from anndata import AnnData


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

    write_spatial_anndata(
        file_path="test_visium.zarr",
        image_axes=["c", "y", "x"],
        image=np.swapaxes(image, 2, 0),
        tables_adata=tables_adata,
        circles_adata=circles_adata,
    )

    return


if __name__ == "__main__":
    write_visium_adata()
