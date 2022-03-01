from abc import ABC
import os
import sqlite3
import stat
import asyncio
import json
from math import floor
from datetime import datetime
from pathlib import Path
from enum import Enum
from cuid import cuid
import logging
from jkolyer.uploader import S3Uploader

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

def dateSinceEpoch(mydate=datetime.now()):
    result = (mydate - datetime(1970, 1, 1)).total_seconds()
    return floor(result)

class UploadStatus(Enum):
    PENDING = 1
    IN_PROGRESS = 2
    COMPLETED = 3
    FAILED = 4

class BaseModel(ABC):
    db_name = 'parallel-file-upload.db'
    db_conn = sqlite3.connect(db_name)
    bucket_name = 'rewotes-pfu-bucket'
    
    @classmethod
    def create_tables(self):
        cursor = self.db_conn.cursor()
        try:
            sqls = FileModel.create_table_sql()
            for sql in sqls: cursor.execute(sql)
            self.db_conn.commit()
            
            sqls = BatchJobModel.create_table_sql()
            for sql in sqls: cursor.execute(sql)
            self.db_conn.commit()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()

    @classmethod
    def run_sql_query(self, sql):
        cursor = self.db_conn.cursor()
        try:
            return cursor.execute(sql).fetchall()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()

    @classmethod
    def run_sql_command(self, sql):
        cursor = self.db_conn.cursor()
        try:
            cursor.execute(sql)
            self.db_conn.commit()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()

    def __init__(self):
        pass

class FileModel(BaseModel):

    @classmethod
    def table_name(cls):
        return 'FileStat'
    
    @classmethod
    def create_table_sql(cls):
        return ["""
        CREATE TABLE IF NOT EXISTS {table_name}
                  ( id TEXT PRIMARY KEY, 
                    created_at INTEGER,
                    file_size INTEGER,
                    last_modified INTEGER,
                    permissions TEXT,
                    file_path TEXT,
                    status INTEGER
                  );
                """.format(
                    table_name=cls.table_name(),
                ),
                f"CREATE UNIQUE INDEX IF NOT EXISTS IdxFilePath ON {cls.table_name()}(file_path)",
                f"CREATE INDEX IF NOT EXISTS IdxStatus ON {cls.table_name()}(status);"]

    @classmethod
    def bootstrap_table(cls):
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
        data = {
            "file_size": self.file_size,
            "last_modified": self.last_modified,
            "permissions": self.permissions,
        }
        return json.dumps(data)

    def _update_status(self, cursor):
        sql = f"UPDATE {self.table_name()} SET status = {self.status} WHERE id = '{self.id}'"
        cursor.execute(sql)
        self.db_conn.commit()

    def start_upload(self, cursor):
        self.status = UploadStatus.IN_PROGRESS.value
        self._update_status(cursor)
        
        completed = self.uploader.upload_file(self.file_path, self.bucket_name, self.id)
        if completed:
            completed = self.uploader.upload_metadata(self.metadata(), self.bucket_name, f"metadata-{self.id}")
        self.upload_complete(cursor) if completed  else self.upload_failed(cursor) 

    def upload_complete(self, cursor):
        self.status = UploadStatus.COMPLETED.value
        self._update_status(cursor)

    def upload_failed(self, cursor):
        self.status = UploadStatus.FAILED.value
        self._update_status(cursor)

    def get_uploaded_file(self):
        return self.uploader.get_uploaded_data(self.bucket_name, self.id)

    def get_uploaded_metadata(self):
        metadata = self.uploader.get_uploaded_data(self.bucket_name, f"metadata-{self.id}")
        return json.loads(metadata)

        
class BatchJobModel(BaseModel):

    def __init__(self, props):
        self.id = props[0]
        self.status = props[1]
        self.created_at = props[2]
        self.root_dir = props[3]

    @classmethod
    def table_name(cls):
        return 'BatchJob'
    
    @classmethod
    def create_table_sql(cls):
        return ["""
        CREATE TABLE IF NOT EXISTS {table_name}
        ( id TEXT PRIMARY KEY, 
        status INTEGER,
        created_at INTEGER,
        root_dir TEXT
        );
        """.format(table_name=cls.table_name()),
                f"CREATE INDEX IF NOT EXISTS IdxCreatedAt ON {cls.table_name()}(created_at);",]

    @classmethod
    def new_record_sql(cls, root_dir):
        return """
        INSERT INTO {table_name}
                  ( id, status, created_at, root_dir )
                  VALUES 
                  ( '{idval}', {status}, {created_at}, '{root_dir}' )
                """.format(
                    table_name=cls.table_name(),
                    idval=cuid(),
                    status=UploadStatus.PENDING.value,
                    created_at=dateSinceEpoch(),
                    root_dir=root_dir,)
    
    @classmethod
    def query_latest(cls):
        sql = f"SELECT * FROM {cls.table_name()} ORDER BY created_at DESC LIMIT 1"
        cursor = cls.db_conn.cursor()
        try:
            result = cursor.execute(sql).fetchall()
            if len(result) == 0: return None
            
            # logger.debug(f"BatchJobModel.query_latest: {result}")
            model = BatchJobModel(result[0])
            return model
        
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()
        return None

    def generate_file_records(self):
        cursor = self.db_conn.cursor()
        file_count = 0
        try:
            for file_path in Path(self.root_dir).rglob('*'):
                fstat = os.stat(file_path)
                fmode = fstat.st_mode
                if stat.S_ISDIR(fmode): continue

                # logger.debug(file_path)
                file_size = fstat.st_size
                last_modified = fstat.st_mtime
                permissions = oct(fstat.st_mode)[-3:]
                status = UploadStatus.PENDING.value

                file_obj = FileModel((
                    cuid(),
                    dateSinceEpoch(),
                    file_size,
                    last_modified,
                    permissions,
                    file_path,
                    status
                ))
                file_obj.save(cursor)
                self.db_conn.commit()
                
                file_count += 1
        
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}")
        finally:
            cursor.close()
        return file_count

    def _fetch_files(self, cursor, page_num, page_size):
        offset = page_num * page_size
        
        # paginate without using sql OFFSET https://gist.github.com/ssokolow/262503
        sql = """
        SELECT * FROM {table_name} 
        WHERE status = {status} AND 
        (id NOT IN ( SELECT id FROM {table_name} ORDER BY file_size ASC LIMIT {offset} ))
        ORDER BY file_size ASC
        LIMIT {page_size}
        """.format(
            table_name = FileModel.table_name(),
            status = UploadStatus.PENDING.value,
            offset = offset,
            page_size = page_size
        )
        results = cursor.execute(sql).fetchall()
        return results
        

    def file_iterator(self, cursor=None):
        _cursor = cursor if cursor else self.db_conn.cursor()

        page_num = 0
        page_size = 10
        try:
            while True:
                results = self._fetch_files(_cursor, page_num, page_size)
                if len(results) == 0: break

                page_num += 1
                for result in results:
                    model = FileModel(result)
                    # logger.debug(f"id = {model.id}")
                    yield model, _cursor
                    
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}")
        finally:
            if cursor is None:
                _cursor.close

    async def async_upload_files(self):
        cursor = self.db_conn.cursor()
        max_concur = 8
        sem = asyncio.Semaphore(max_concur)

        async def task_wrapper(model, cursor):
            # logger.debug(f"task_wrapper:  model = {model.file_path}")
            try:
                model.start_upload(cursor)
            finally:
                sem.release()
                
        for model, cursor in self.file_iterator(cursor):
            await sem.acquire()
            asyncio.create_task(task_wrapper(model, cursor))

        # wait for all tasks to complete
        for i in range(max_concur):
            await sem.acquire()
        cursor.close()

    def reset_file_status(self):
        cursor = self.db_conn.cursor()
        try:
            sql = f"UPDATE {FileModel.table_name()} SET status = {UploadStatus.PENDING.value}"
            cursor.execute(sql)
            self.db_conn.commit()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; {sql}")
        finally:
            cursor.close
        
        
        
