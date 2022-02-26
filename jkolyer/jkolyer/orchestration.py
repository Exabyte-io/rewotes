import sqlite3
import logging
from jkolyer.models import FileModel, JobModel

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class Orchestration:
    def __init__(self):
        pass

    def disconnect_db(self):
        if self.db_conn is None: return
        self.db_conn.close()
        self.db_conn = None
        logger.info("The SQLite connection is closed")
    
    def connect_db(self):
        try:
            self.db_conn = sqlite3.connect('parallel-file-upload.db')
            cursor = self.db_conn.cursor()

            sqlite_select_Query = "select sqlite_version();"
            cursor.execute(sqlite_select_Query)
            record = cursor.fetchall()
            cursor.close()

        except sqlite3.Error as error:
            self.db_conn = None
            logger.error(f"Error while connecting to sqlite: {error}")
            
    def create_tables(self):
        if self.db_conn is None: return
        cursor = self.db_conn.cursor()
        try:
            sqls = FileModel.create_table_sql()
            for sql in sqls: cursor.execute(sql)
            self.db_conn.commit()
            
            sqls = JobModel.create_table_sql()
            for sql in sqls: cursor.execute(sql)
            self.db_conn.commit()

        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")

        finally:
            cursor.close()

    def run_sql(self, sql):
        if self.db_conn is None: return
        cursor = self.db_conn.cursor()
        try:
            return cursor.execute(sql).fetchall()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()
        

