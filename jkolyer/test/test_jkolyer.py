"""
Tests for `jkolyer` module.
"""
import pytest
from jkolyer.orchestration import Orchestration
from jkolyer.models import FileModel, JobModel


class TestJkolyer(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    
@pytest.fixture
def orchestration():
    orch = Orchestration()
    orch.connect_db()
    return orch
    

class TestOrchestration(TestJkolyer):
    def test_create(self, orchestration):
        assert orchestration != None

    def test_connectDb(self, orchestration):
        orchestration.connect_db()
        assert orchestration.db_conn is not None
        orchestration.disconnect_db()
        assert orchestration.db_conn is None

    def test_create_tables(self, orchestration):
        orchestration.create_tables()
        sql = f"SELECT name FROM sqlite_master WHERE type='table' AND name='{FileModel.table_name()}'"
        result = orchestration.run_sql(sql)
        assert result[0][0] == FileModel.table_name()
        sql = f"SELECT name FROM sqlite_master WHERE type='table' AND name='{JobModel.table_name()}'"
        result = orchestration.run_sql(sql)
        assert result[0][0] == JobModel.table_name()
        
