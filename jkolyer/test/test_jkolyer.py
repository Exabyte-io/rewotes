import os
import pytest
import sqlite3
import boto3
from moto import mock_s3
from threading import Thread
import logging

for name in logging.Logger.manager.loggerDict.keys():
    if ('boto' in name) or \
       ('urllib3' in name) or \
       ('boto3' in name) or \
       ('botocore' in name) or \
       ('nose' in name):
        logging.getLogger(name).setLevel(logging.CRITICAL)
logging.getLogger('s3transfer').setLevel(logging.CRITICAL)                    

from jkolyer.models.base_model import BaseModel, UploadStatus
from jkolyer.models.batch_model import BatchJobModel, parallel_upload_files
from jkolyer.models.file_model import FileModel
from jkolyer.uploader import S3Uploader

    
@pytest.fixture
def batch_job():
    return BatchJobModel.new_instance('./samples')
    
class TestJkolyer(object):
    @classmethod
    def setup_class(cls):
        S3Uploader.set_boto3_client(S3Uploader.s3_mock())

    @classmethod
    def teardown_class(cls):
        pass

class TestTables(TestJkolyer):
    def test_create_tables(self):
        FileModel.create_tables()
        BatchJobModel.create_tables()
        sql = f"SELECT name FROM sqlite_master WHERE type='table' AND name='{FileModel.table_name()}'"
        result = BaseModel.run_sql_query(sql)
        assert result[0][0] == FileModel.table_name()

class TestBatchJob(TestJkolyer):
    
    def test_create_table(self):
        sql = f"SELECT name FROM sqlite_master WHERE type='table' AND name='{BatchJobModel.table_name()}'"
        result = BaseModel.run_sql_query(sql)
        assert result[0][0] == BatchJobModel.table_name()

    def test_create_batch(self, batch_job):
        result = BatchJobModel.query_latest()
        assert result is not None

class TestFileModel(TestJkolyer):
    
    @classmethod
    def setup_class(cls):
        TestJkolyer.setup_class()
        FileModel.bootstrap_table()

    def test_create_file_records(self):
        batch = BatchJobModel.query_latest()
        file_count = batch.generate_file_records()
        
        cursor = BaseModel.db_conn.cursor()
        result = cursor.execute(f"SELECT COUNT(*) FROM {FileModel.table_name()}").fetchall()
        assert result[0][0] == file_count

        # ensure no duplicates are created
        new_file_count = batch.generate_file_records()
        result = cursor.execute(f"SELECT COUNT(*) FROM {FileModel.table_name()}").fetchall()
        assert result[0][0] == file_count
        
        cursor.close()

    def test_file_upload(self):
        model = FileModel.fetch_record(UploadStatus.PENDING.value)
        assert model is not None
        cursor = FileModel.db_conn.cursor()
        model.start_upload(cursor)
        assert model.status == UploadStatus.COMPLETED.value
        model2 = FileModel.fetch_record(UploadStatus.COMPLETED.value)
        assert model2 is not None
        assert model2.id == model.id
        cursor.close()
        
        file_contents = model.get_uploaded_file()
        assert file_contents is not None
        metadata = model.get_uploaded_metadata()
        assert metadata is not None
        
    def test_batch_uploads_sequential(self):
        batch = BatchJobModel.query_latest()
        
        for file_model, cursor in batch.file_iterator():
            assert file_model.status == UploadStatus.PENDING.value
            file_model.start_upload(cursor)
            assert file_model.status == UploadStatus.COMPLETED.value
        
class TestAsyncFileModel(TestJkolyer):
    
    @classmethod
    def setup_class(cls):
        batch = BatchJobModel.query_latest()
        # reset the file records
        batch.reset_file_status()

    @pytest.mark.asyncio            
    async def test_batch_uploads_async(self):
        batch = BatchJobModel.query_latest()
        await batch.async_upload_files()

        for file_model, cursor in batch.file_iterator():
            assert file_model.status == UploadStatus.COMPLETED.value

class TestParallelFileModel(TestJkolyer):
    
    @classmethod
    def setup_class(cls):
        batch = BatchJobModel.query_latest()
        batch.reset_file_status()

    def test_batch_uploads_parallel(self): 
        batch = BatchJobModel.query_latest()
        parallel_upload_files(batch)
        
        for file_model, cursor in batch.file_iterator():
            assert file_model.status == UploadStatus.COMPLETED.value
        
