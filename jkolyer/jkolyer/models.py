import os
import sqlite3
import stat
from math import floor
from datetime import datetime
from pathlib import Path
from enum import Enum
from cuid import cuid
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

def dateSinceEpoch(mydate=datetime.now()):
    result = (mydate - datetime(1970, 1, 1)).total_seconds()
    return floor(result)

class BaseModel:
    db_name = 'parallel-file-upload.db'
    db_conn = sqlite3.connect(db_name)
    
    @classmethod
    def create_tables(self):
        cursor = self.db_conn.cursor()
        try:
            sqls = FileModel.create_table_sql()
            for sql in sqls: cursor.execute(sql)
            self.db_conn.commit()
            
            sqls = UploadJobModel.create_table_sql()
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
                    file_path TEXT
                  );
                """.format(
                    table_name=cls.table_name(),
                ),
                f"CREATE UNIQUE INDEX IF NOT EXISTS IdxFilePath ON {cls.table_name()}(file_path)",
                ]

    def __init__(self, file_size, last_modified, permissions, file_path):
        self.id = cuid()
        self.created_at = dateSinceEpoch()
        self.file_size = file_size
        self.last_modified = last_modified
        self.permissions = permissions
        self.file_path = file_path

    def save(self, cursor):
        sql = """
            INSERT OR IGNORE INTO {table_name}
                  ( id, created_at, file_size, last_modified, permissions, file_path )
                  VALUES 
                  ( '{id}', {created_at}, {file_size}, {last_modified}, '{permissions}', '{file_path}' )
                """.format(
                    table_name=self.__class__.table_name(),
                    id=self.id,
                    created_at=self.created_at,
                    file_size=self.file_size,
                    last_modified=self.last_modified,
                    permissions=self.permissions,
                    file_path=self.file_path
                )
        cursor.execute(sql)
        
class UploadJobModel(BaseModel):
    @classmethod
    def table_name(cls):
        return 'UploadJob'
    
    @classmethod
    def create_table_sql(cls):
        return ["""
        CREATE TABLE IF NOT EXISTS {table_name}
                  ( id TEXT PRIMARY KEY, 
                    batch_id TEXT,
                    file_id TEXT,
                    status INTEGER,
                    created_at INTEGER,
                    FOREIGN KEY (file_id) REFERENCES {file_table_name}(id),
                    FOREIGN KEY (batch_id) REFERENCES {batch_table_name}(id)
                  );
               """.format(
                   table_name=cls.table_name(),
                   file_table_name=FileModel.table_name(),
                   batch_table_name=BatchJobModel.table_name()
               ),
                f"CREATE INDEX IF NOT EXISTS IdxJobFile ON {cls.table_name()}(file_id);",
                f"CREATE INDEX IF NOT EXISTS IdxBatch ON {cls.table_name()}(batch_id);",
                f"CREATE INDEX IF NOT EXISTS IdxStatus ON {cls.table_name()}(status);"]

class BatchStatus(Enum):
    PENDING = 1
    IN_PROGRESS = 2
    COMPLETED = 3
    FAILED = 4

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
                f"CREATE INDEX IF NOT EXISTS IdxCreatedAt ON {cls.table_name()}(created_at);",
                ]

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
                    status=BatchStatus.PENDING.value,
                    created_at=dateSinceEpoch(),
                    root_dir=root_dir,)
    
    @classmethod
    def query_latest(cls):
        sql = f"SELECT * FROM {cls.table_name()} ORDER BY created_at DESC LIMIT 1"
        cursor = cls.db_conn.cursor()
        try:
            result = cursor.execute(sql).fetchall()
            if len(result) == 0: return None
            
            logger.debug(f"BatchJobModel.query_latest: {result}")
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

                logger.debug(file_path)
                file_size = fstat.st_size
                last_modified = fstat.st_mtime
                permissions = stat.S_IMODE(fmode)

                file_obj = FileModel(file_size, last_modified, permissions, file_path)
                file_obj.save(cursor)
                self.db_conn.commit()
                
                file_count += 1
        
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}")
        finally:
            cursor.close()
        return file_count
    
