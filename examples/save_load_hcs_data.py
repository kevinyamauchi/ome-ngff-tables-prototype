import napari
import os
import shutil
import squidpy as sq
from ngff_tables_prototype.writer import write_spatial_anndata
from ngff_tables_prototype.reader import load_to_napari_viewer
import numpy as np
import pandas as pd
from pathlib import Path
import anndata as ad

# Instructions
# Download files from here:
# https://drive.google.com/drive/folders/1TWMcKSQuq_ZEkTPPCLOO4sE6oMlhSGu_?usp=sharing
# Place it in the ome-ngff-tables-prototype folder and run the script from
# the example folder
# Unzip the OME-Zarr file
output_fpath = Path("../HCS_Example_Data/UZH_2x2_HCS.ome.zarr")
feature_path = Path("../HCS_Example_Data/FeatureData")

def write_hcs_table() -> None:
    # Loop over the 4 sites
    # Check that data was downloaded
    assert len(os.listdir(feature_path)) == 4, 'Feature data should contain \
                                            for files. Check instructions above'
    assert os.path.isdir(output_fpath), 'OME Zarr folder does not exist. Was it unzipped?'

    # TODO: Check if feature data has already been written to the OME-Zarr file

    # Loop through all 4 fields of view,
    # load the feature data & save to OME-Zarr file
    for site in range(4):
        df = pd.read_csv(feature_path / f'HCS_FeatureData_Site{site}.csv')
        df_x = df.iloc[:, :-2]
        features_ad = ad.AnnData(df_x.loc[:, df_x.columns != 'Label'])
        features_ad.obs = pd.concat([df['Label'], df.iloc[:, -2:]], axis = 1)
        write_spatial_anndata(
            file_path=output_fpath / 'B' / '03' / str(site),
            tables_adata=features_ad,
            tables_region="labels/label_image",
            tables_instance_key="Label",
        )

    return


if __name__ == "__main__":
    write_hcs_table()
    print(output_fpath)
    viewer = load_to_napari_viewer(
        file_path=output_fpath / 'B' / '03' / '0',
        groups=["labels/label_image", "tables/regions_table"],
    )
    napari.run()
