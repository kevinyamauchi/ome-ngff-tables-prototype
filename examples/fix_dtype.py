
import sys
import zarr
import numpy as np
import argparse

# convert zarr dypes from int64 (i8) to int32 (i4) to allow data
# to be loaded via JavaScript using zarr.js

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='Path to e.g. "test_segment.zarr/tables/regions_table/obs"')
    args = parser.parse_args(argv)
    zarr_path = args.path


    name_list = []
    data_list = []
    chunk_list = []
    # store = zarr.DirectoryStore(zarr_path, dimension_separator="/")
    with zarr.open(zarr_path, "a") as f:
        col_names = []
        f.visit(lambda x: col_names.append(x))
        print('col_names', col_names)
        for ds_name in col_names:
            data = f[ds_name]
            # if we have an int64 array...
            if hasattr(data, "dtype") and data.dtype == np.int64:
                print("fix dtype...", ds_name)
                data = f[ds_name][:].astype("int32")
                chunks = f[ds_name].chunks
                # delete and replace...
                del f[ds_name]
                name_list.append(ds_name)
                data_list.append(data)
                chunk_list.append(chunks)
        
        for ii, ds_name in enumerate(name_list):
            print("write ii, ds_name", ii, ds_name)
            data = data_list[ii]
            chunks = chunk_list[ii]
            f.create_dataset(ds_name, data=data, chunks=chunks)


if __name__ == "__main__":
    main(sys.argv[1:])
