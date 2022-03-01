"""FileModel: the data model for the FileStat database table.

Provides SQL wrapper around file metadata and upload status.
"""
import sqlite3
import json
import logging
from jkolyer.uploader import S3Uploader
from jkolyer.models.base_model import BaseModel, UploadStatus

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class FileModel(BaseModel):

    @classmethod
    def table_name(cls):
        """Returns the SQL table name 'FileStat'
        :return: string 
        """
        return 'FileStat'
    
    @classmethod
    def create_table_sql(cls):
        """All the sql create scripts needed by file objects 
           for tables and indices.  
           Does nothing if the tables/indices already exist.
        :return: string[] SQL statements
        """
        sql = """CREATE TABLE IF NOT EXISTS {table_name}
           ( id TEXT PRIMARY KEY, 
            created_at INTEGER,
            file_size INTEGER,
            last_modified INTEGER,
            permissions TEXT,
            file_path TEXT,
            status INTEGER
            );""".format(table_name=cls.table_name())
        return [sql,
                f"CREATE UNIQUE INDEX IF NOT EXISTS IdxFilePath ON \
                   {cls.table_name()}(file_path)",
                f"CREATE INDEX IF NOT EXISTS IdxStatus ON \
                   {cls.table_name()}(status);"]

    @classmethod
    def bootstrap_table(cls):
        """Drops and recreates SQL tables
        :return: None
        """
        cursor = cls.db_conn.cursor()
        try:
            cursor.execute(f"DROP TABLE IF EXISTS {cls.table_name()}")
            cls.db_conn.commit()
            for sql in cls.create_table_sql(): cursor.execute(sql)
            cls.db_conn.commit()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()

    @classmethod
    def fetch_record(cls, status):
        """Retrieves the most recently-created instance with the given status.
        :param status (`UploadStatusEnum`): target status
        :return: FileModel instance if found, otherwise None
        """
        sql = f"SELECT * FROM {FileModel.table_name()} WHERE status = {status} ORDER BY created_at DESC LIMIT 1"
        cursor = cls.db_conn.cursor()
        try:
            result = cursor.execute(sql).fetchone()
            return FileModel(result) if result is not None else None
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()
        return None

    def __init__(self, *args):
        """Instance constructor, setting table properties, and local `S3Uploader` instance.
        :param args: tuple of values ordered as in create table script
        """
        tpl = args[0]
        self.id = tpl[0]
        self.created_at = tpl[1]
        self.file_size = tpl[2]
        self.last_modified = tpl[3]
        self.permissions = tpl[4]
        self.file_path = tpl[5]
        self.status = tpl[6]
        self.uploader = S3Uploader()

    def save(self, cursor):
        """Saves the receiver's properties into the database using INSERT OR IGNORE statement.
           Will throw exception on error.
        :param cursor: active cursor to execute SQL
        :return: None
        """
        sql = """
            INSERT OR IGNORE INTO {table_name}
                  ( id, created_at, file_size, last_modified, permissions, file_path, status )
                  VALUES 
                  ( '{id}', {created_at}, {file_size}, {last_modified}, '{permissions}', '{file_path}', {status} )
                """.format(
                    table_name=self.__class__.table_name(),
                    id=self.id,
                    created_at=self.created_at,
                    file_size=self.file_size,
                    last_modified=self.last_modified,
                    permissions=self.permissions,
                    file_path=self.file_path,
                    status=self.status
                )
        cursor.execute(sql)

    def metadata(self):
        """Data structure used to store file metadata in storage provider.
           Properties include `file_size`, `last_modified`, and `permissions`.
           Rendered as string for storage purposes.
        :return: string JSON-formatted using `json.dumps`
        """
        data = {
            "file_size": self.file_size,
            "last_modified": self.last_modified,
            "permissions": self.permissions,
        }
        return json.dumps(data)

    def _update_status(self, cursor):
        """Convenience method for SQL UPDATE of the `status` property.
           Executes SQL and commits.  Throws exception on error.
        :param cursor: used for SQL execution
        :return: None
        """
        sql = f"UPDATE {self.table_name()} SET status = {self.status} WHERE id = '{self.id}'"
        cursor.execute(sql)
        self.db_conn.commit()

    def start_upload(self, cursor):
        """Status state change to initiate upload.  Calls `_update_status`.
           Invokes `upload_file` and `upload_metadata` on the uploader property.
           If either upload fails, calls `upload_failed`; otherwise calls `upload_complete`
        :param cursor: used for SQL execution
        :return: None
        """
        self.status = UploadStatus.IN_PROGRESS.value
        self._update_status(cursor)
        
        completed = self.uploader.upload_file(self.file_path, self.bucket_name, self.id)
        if completed:
            completed = self.uploader.upload_metadata(self.metadata(), self.bucket_name, f"metadata-{self.id}")
        self.upload_complete(cursor) if completed  else self.upload_failed(cursor) 

    def upload_complete(self, cursor):
        """Status state change to success upload completion.  Calls `_update_status`.
        :param cursor: used for SQL execution
        :return: None
        """
        self.status = UploadStatus.COMPLETED.value
        self._update_status(cursor)

    def upload_failed(self, cursor):
        """Status state change to failed upload.  Calls `_update_status`.
        :param cursor: used for SQL execution
        :return: None
        """
        self.status = UploadStatus.FAILED.value
        self._update_status(cursor)

    def get_uploaded_file(self):
        """For testing purposes, fetches the uploaded file from object storage
        :return: binary string: the uploaded file bytes or None
        """
        return self.uploader.get_uploaded_data(self.bucket_name, self.id)

    def get_uploaded_metadata(self):
        """For testing purposes, fetches the uploaded metadata from object storage
        :return: dict: the uploaded metadata or None
        """
        metadata = self.uploader.get_uploaded_data(self.bucket_name, f"metadata-{self.id}")
        return json.loads(metadata)

    def parallel_dto_string(self):
        """For parallel uploading with multiprocessing module,
           provides data transfer object needed for uploading
           in a separate process: `id`, `file_path`, `metadata`,
           `bucket_name`, `status`.
        :return: string: the JSON string of properties
        """
        dto = {
            "id": self.id,
            "file_path": self.file_path,
            "metadata": self.metadata(),
            "bucket_name": self.bucket_name,
            "status": self.status,
        }
        return json.dumps(dto)


