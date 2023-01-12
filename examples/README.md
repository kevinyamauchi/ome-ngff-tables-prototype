
To view sample data in Vitessce, generate the data:

    $ cd examples
    $ python save_load_squidpy_segment.py


Then need to convert test_segment.zarr/.zattrs to a v0.3 image and add 'omero' metadata for channels:

.zattrs:

```
{
    "multiscales": [
        {
            "axes": [
                "c", "y", "x"
            ],
            "datasets": [
                {"path": "0"},
                {"path": "1"},
                {"path": "2"},
                {"path": "3"},
                {"path": "4"}
            ],
            "name": "/",
            "version": "0.3"
        }
    ],
    "omero": {
        "id": 1,                          
        "name": "example",               
        "version": "0.3",                
        "channels": [
            {
                "active": true,
                "coefficient": 1,
                "color": "FF0000",
                "family": "linear",
                "inverted": false,
                "label": "LaminB1",
                "window": {
                    "end": 255,
                    "max": 255,
                    "min": 0,
                    "start": 0
                }
            },
            {
                "active": true,
                "coefficient": 1,
                "color": "00FF00",
                "family": "linear",
                "inverted": false,
                "label": "Green",
                "window": {
                    "end": 255,
                    "max": 255,
                    "min": 0,
                    "start": 0
                }
            },
            {
                "active": true,
                "coefficient": 1,
                "color": "0000FF",
                "family": "linear",
                "inverted": false,
                "label": "Blue",
                "window": {
                    "end": 255,
                    "max": 255,
                    "min": 0,
                    "start": 0
                }
            }
        ],
        "rdefs": {
            "defaultT": 0,
            "defaultZ": 0,
            "model": "color"
        }
    }
}
```

Serve the data from `localhost:9000`, e.g. using https://github.com/http-party/http-server:

    $ http-server ./ --cors -p 9000

In Firefox (not Chrome)...(https://github.com/vitessce/vitessce/issues/1089)
go to http://vitessce.io/#?edit=false&url=http://localhost:9000/vitessce_anndata_config.json
