import pytest
import sqlite3
from jkolyer.orchestration import Orchestration
from jkolyer.models import BaseModel, BatchJobModel, FileModel, UploadJobModel


class TestJkolyer(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    
@pytest.fixture
def batch_job():
    sql = BatchJobModel.new_record_sql('./samples')
    BaseModel.run_sql_command(sql)
    

class TestTables(TestJkolyer):
    def test_create_tables(self):
        BaseModel.create_tables()
        sql = f"SELECT name FROM sqlite_master WHERE type='table' AND name='{FileModel.table_name()}'"
        result = BaseModel.run_sql_query(sql)
        assert result[0][0] == FileModel.table_name()
        
        sql = f"SELECT name FROM sqlite_master WHERE type='table' AND name='{UploadJobModel.table_name()}'"
        result = BaseModel.run_sql_query(sql)
        assert result[0][0] == UploadJobModel.table_name()

        
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
        cursor = BaseModel.db_conn.cursor()
        cursor.execute(f"DROP TABLE IF EXISTS {FileModel.table_name()}")
        BaseModel.db_conn.commit()
        for sql in FileModel.create_table_sql(): cursor.execute(sql)
        BaseModel.db_conn.commit()
        cursor.close()

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
        
