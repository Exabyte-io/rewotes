# Multi-threaded file uploader (Backend)

Can be run using `./run` in the directory meant for syncing.
Can be tested with `./test`.

AWS Credentials can be stored in the environment variables
`AWS_ID`, `AWS_KEY`, and `AWS_BUCKET`, to specify the location of the uploads.


The program has the following files.

- main:
    This is responsible for scheduling the sync every hour.

- scanner:
    Responsible for reading all the files in the given directory.
    Fills the buffer up to certain chunk size and uploads the batch.

- buffer:
    Responsible for checking which files have been updated.
    Updates the metadata.

- uploader:
    Responsible for uploading given files to cloud storage.

- file:
    Wrapper for Python file object.
    Stores the hash, size and filename as class variables.

- client:
    Responsible for direct interaction with cloud storage.
    Handles credentials and uploading files.


__Notes on improvements.__
The credential management system and user permissions could be improved.

The uploader does not handle the situation of single files that are very large. This could be improved by allowing some sort of chunking of the file, and uploading each chunk individually.

Adding more unit-tests would be important.

The Client module should support multiple cloud storage types. This could be done by changing it to an abstract class that needs the `metadata()`, `save_metadata()`, and `upload()` methods. It could then be implemented by different versions (Google Drive, Azue, AWS S3) for example.
