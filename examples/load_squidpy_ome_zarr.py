import napari
from ngff_tables_prototype.reader import load_to_napari_layer_data

file_path = "../notebooks/test_anndata.zarr"
table_group = "table"

layer_data = load_to_napari_layer_data(file_path=file_path, table_group_key=table_group)

viewer = napari.Viewer()
for layer in layer_data:
    viewer._add_layer_from_data(*layer)

napari.run()

