# Adding New Features to the Model

To add a new feature to the model, use the following procedure.

Following this procedure will take approximately 5 to 10 minutes.

### Step 1: Modify the Material Protocol Buffer

#### 1.1

The `Material` class is a thin wrapper around a protocol buffer defined in `pkg/material/material.proto`.
Add a new field to the `Material` proto definition in this file to be able to store your desired new feature.
For a syntax guide on protocol buffers see: [https://developers.google.com/protocol-buffers/docs/proto3](https://developers.google.com/protocol-buffers/docs/proto3)

#### 1.2

Compile the `material.proto` file to `material_pb2.py` in the same directory using the following.

```bash
cd pkg/material
protoc -I=. --python_out=. material.proto
```

The `protoc` protocol buffer compiler can be downloaded from [https://github.com/protocolbuffers/protobuf/releases/latest](https://github.com/protocolbuffers/protobuf/releases/latest).

Instructions for using the `protoc` protocol buffer compiler can be found at [https://developers.google.com/protocol-buffers/docs/pythontutorial#compiling-your-protocol-buffers](https://developers.google.com/protocol-buffers/docs/pythontutorial#compiling-your-protocol-buffers).

#### Conclusion

The new field will automatically be available as a data member in the `Material` Python class defined in `pkg/material/Material.py`.

### Step 2: Modify the download method

The `Downloader.download` method, defined in `pkg/downloader/Downloader.py`, specifies which fields to download from `materialsproject.org`.
Add the name of your field, as it appears in the Materials Project API, to the list of fields to be downloaded. As of the time of writing, the Materials Project API does not appear to be fully documented, so determining the exact name of your field is left as an exercise for the reader. The field names appear to be lowercase English words separated by underscores, although there is no guarantee that all field names follow this pattern.

Later in the `Downloader.download` method, add an assignment from the data structure received from the Materials Project API to the `Material` class.

Examples for both of these tasks exist in the method's source code already.
When in doubt, copy/paste.

### Step 3: Modify the to_numpy_double_array method

The `Material.to_numpy_double_array` method, defined in `pkg/material/Material.py`, defines the encoding of a `Material` class into an array of numbers.
At the end of this function, just before the return statement, read from `self.<your_protocol_buffer_field_name>`, then append a `float` or `float`s
based on the read value.

Examples for this task exist in the method's source code.
When in doubt, copy/paste.

### Conclusion

At this point, the field will be fully integrated into the package.
The field can be read from and written to the `Material` class,
the field can be downloaded from the Materials Project,
and the field will be incorporated when generating training, testing, and prediction rows
for the machine learning model.
