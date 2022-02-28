import sqlite3
from jkolyer.models import FileModel, UploadJobModel, BatchJobModel
import logging
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
            self.db_conn = sqlite3.connect(BatchJobModel.db_name())
            cursor = self.db_conn.cursor()

            sqlite_select_Query = "select sqlite_version();"
            cursor.execute(sqlite_select_Query)
            record = cursor.fetchall()
        except sqlite3.Error as error:
            self.db_conn = None
            logger.error(f"Error while connecting to sqlite: {error}")
        finally:
            cursor.close()
