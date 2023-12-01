"""BatchJobModel: the data model for the BatchJob database table.

Provides SQL wrapper around upload status for a set of files.
"""
import os
import stat
from pathlib import Path
import sqlite3
import asyncio
import json
from multiprocessing import cpu_count, Process, Queue, Semaphore, Value
from cuid import cuid
import logging
import time

from jkolyer.models.base_model import BaseModel, UploadStatus, dateSinceEpoch
from jkolyer.models.file_model import FileModel
from jkolyer.uploader import S3Uploader

for name in logging.Logger.manager.loggerDict.keys():
    if ('boto' in name) or \
       ('urllib3' in name) or \
       ('boto3' in name) or \
       ('botocore' in name) or \
       ('nose' in name):
        logging.getLogger(name).setLevel(logging.CRITICAL)
logging.getLogger('s3transfer').setLevel(logging.CRITICAL)                    

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)
        
class BatchJobModel(BaseModel):

    @classmethod
    def table_name(cls):
        """Returns the SQL table name 'BatchJob'
        :return: string 
        """
        return 'BatchJob'
    
    @classmethod
    def create_table_sql(cls):
        """All the sql create scripts needed by file objects 
           for tables and indices.  
           Does nothing if the tables/indices already exist.
        :return: string[] SQL statements
        """
        return [
            """CREATE TABLE IF NOT EXISTS {table_name}
                ( id TEXT PRIMARY KEY, 
                status INTEGER,
                created_at INTEGER,
                root_dir TEXT );
            """.format(table_name=cls.table_name()),
            f"CREATE INDEX IF NOT EXISTS IdxCreatedAt ON {cls.table_name()}(created_at);",
        ]

    @classmethod
    def new_record_sql(cls, root_dir):
        """Provides the SQL for creating a new record with default values.
        :param root_dir: the file directories root path
        :return: string SQL INSERT statement
        """
        return """INSERT INTO {table_name}
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
        """Fetches the most recent record from the database
        :return: BatchModel: latest instance or None
        """
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

    @classmethod
    def new_instance(cls, root_dir):
        """Creates and return new instance
        :param root_dir: the file directories root path
        :return: BatchModel: the instance
        """
        sql = cls.new_record_sql(root_dir)
        cls.run_sql_command(sql)
        return cls.query_latest()
    
    def __init__(self, props):
        """Instance constructor, setting table properties
        :param args: tuple of values ordered as in create table script
        """
        self.id = props[0]
        self.status = props[1]
        self.created_at = props[2]
        self.root_dir = props[3]

    def generate_file_records(self):
        """Loads all files from receiver's `root_dir` and 
           creates a `FileModel` instance for each,
           which then saves a new database record.
        :return: int: the count of file records created
        """
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
        """Convenience method for retrieving a page of database
           records for `FileModel`.  Orders by `file_size` (ascending).
        :param cursor: used for SQL execution
        :param page_num: page number for query offset
        :param page_size: records per page
        :return: tuple[]: an array of tuples for FileModel
        """
        offset = page_num * page_size
        
        # paginate without using sql OFFSET https://gist.github.com/ssokolow/262503
        sql = """
        SELECT * FROM {table_name} 
        WHERE status = {status} AND 
        (id NOT IN ( 
           SELECT id FROM {table_name} 
           ORDER BY file_size ASC LIMIT {offset}
        ))
        ORDER BY file_size ASC
        LIMIT {page_size}
        """.format(
            table_name = FileModel.table_name(),
            status = UploadStatus.PENDING.value,
            offset = offset,
            page_size = page_size
        )
        return cursor.execute(sql).fetchall()

    def file_iterator(self, cursor=None):
        """Generator method for iterating over a page of `FileModel` data.
           Yields to caller with model instance and cursor.  Page size is 10 records.
        :param cursor: for SQL execution; creates cursor if not provided
        :return: None
        """
        _cursor = cursor if cursor else self.db_conn.cursor()

        page_num = 0
        page_size = 10
        try:
            while True:
                results = self._fetch_files(_cursor, page_num, page_size)
                if len(results) == 0: break

                page_num += 1
                for result in results:
                    yield FileModel(result), _cursor
                    
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}")
        finally:
            """only close the cursor if created in this method"""
            if cursor is None:
                _cursor.close

    def reset_file_status(self):
        """Resets all the `FileModel` status values to `UploadStatus.PENDING`.
           Useful for testing or restarting a previous batch upload.
        :return: None
        """
        cursor = self.db_conn.cursor()
        try:
            sql = f"UPDATE {FileModel.table_name()} SET status = {UploadStatus.PENDING.value}"
            cursor.execute(sql)
            self.db_conn.commit()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; {sql}")
        finally:
            cursor.close
        
    async def async_upload_files(self):
        """Performs asynchronous file upload across all pending `FileModel` instances.
           Uses `asyncio` module.  Maximum 8 concurrent jobs.  
        :return: None
        """
        cursor = self.db_conn.cursor()
        max_concur = 8
        sem = asyncio.Semaphore(max_concur)

        async def task_wrapper(model, cursor):
            try:
                model.start_upload(cursor)
            finally:
                sem.release()
                
        for file_model, _cursor in self.file_iterator(cursor):
            """ this blocks until concurrent jobs drops below max """
            await sem.acquire()
            asyncio.create_task(task_wrapper(file_model, _cursor))

        # wait for all tasks to complete
        for i in range(max_concur):
            await sem.acquire()
        cursor.close()

        
MAX_NICENESS = 19

def parallel_upload(file_dto_string, queue, sema):
    """Perform upload command in multiprocessing mode.  This runs 
       seperate from the main process.  Sets process nice level to 19.
       Sets the uploader endpoint based on dto.
    :param file_dto_string:  JSON-formatted string as FileModel dto
    :param queue: inter-process queue to share upload results
    :param sema: semaphora to limit number of active processes
    :return: None
    """
    os.nice(MAX_NICENESS)
    file_dto = json.loads(file_dto_string)

    # logger.debug(f"parallel_upload:  file_dto_string = {file_dto_string}")
    S3Uploader.set_endpoint_url(file_dto["endpoint_url"])
    
    uploader = S3Uploader()
    
    metadata = file_dto["metadata"]
    
    completed = uploader.upload_file(
        file_dto["file_path"],
        file_dto["bucket_name"],
        file_dto["id"],
        metadata["file_size"],
    )
    if completed:
        completed = uploader.upload_metadata(
            json.dumps(metadata), file_dto["bucket_name"], f"metadata-{file_dto['id']}"
        )
        
    file_dto["status"] = UploadStatus.COMPLETED.value if completed else UploadStatus.FAILED.value
    file_dto_string = json.dumps(file_dto)
    
    queue.put(file_dto_string)
    
    # `release` will add 1 to `sema`, allowing other 
    # processes blocked on it to continue
    sema.release()

    
def parallel_upload_files(batch_model):
    """Uploads files using multiprocessing.  
       Limits concurrent process count to half `multiprocessing.cpu_count`.
    :param batch_model: the BatchJobModel instance
    :return: None
    """
    concurrency = cpu_count() // 2
    """ used to limit the total number of processes spawned"""
    sema = Semaphore(concurrency)
    active_processes = {}
    queue = Queue()
    """temporary store of file_models in progress"""
    file_models_progress = {}

    def _read_queue():
        """Handle data in the queue which is 
           populated by the `parallel_upload` function
           run in separate process.  Here we update database
           using `FileModel` state change methods.
        """
        cursor = BatchJobModel.db_conn.cursor()
        try: 
            while not queue.empty():
                dto = queue.get()
                dto = json.loads(dto)
                fmodel = file_models_progress[dto["id"]]
                try:
                    if dto["status"] == UploadStatus.COMPLETED.value:
                        fmodel.upload_complete(cursor)
                    elif dto["status"] == UploadStatus.FAILED.value:
                        fmodel.upload_failed(cursor)
                        
                    del file_models_progress[fmodel.id]
                    del active_processes[fmodel.id]
                    
                except sqlite3.Error as error:
                    logger.error(f"Error running sql: {error}")
        finally:
            cursor.close()

    for file_model, _cursor in batch_model.file_iterator():
        """once max processes are running, the following `acquire` call
           will block the main process since `sema` has been reduced to 0
        """
        sema.acquire()

        """ cache of `FileModel` used when queue is processed """
        file_models_progress[file_model.id] = file_model
        
        dtoStr = file_model.parallel_dto_string()
        proc = Process(target=parallel_upload, args=(dtoStr, queue, sema))
        active_processes[file_model.id] = proc
        proc.start()
        """ check if queue has any data (also clears active_processes and file_models_progress) """
        _read_queue()

    # inside main process, wait for all processes to finish
    for proc in active_processes.values():
        proc.join()
        
    _read_queue()

