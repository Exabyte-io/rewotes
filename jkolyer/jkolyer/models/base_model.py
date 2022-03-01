""" Abstract superclass for SQL data wrappers.  Includes
    utility functions and classes.
"""
from abc import ABC, abstractmethod
import sqlite3
from datetime import datetime
from enum import Enum
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

def dateSinceEpoch(mydate=datetime.now()):
    """Seconds between given date and epoch date 1970-01-01
    :param mydate: source date or now
    :return: float: seconds since epoch date
    """
    result = (mydate - datetime(1970, 1, 1)).total_seconds()
    return result

class UploadStatus(Enum):
    """Enumerator for upload status state changes.
    """
    FAILED = -1
    PENDING = 0
    IN_PROGRESS = 1
    COMPLETED = 2

class BaseModel(ABC):
    """Holds class properties for database and object storage
    :param db_name: SQLite database name
    :param db_conn: database connection provided by `sqlite3`
    :param bucket_name: object storage bucket name
    """
    db_name = 'parallel-file-upload.db'
    db_conn = sqlite3.connect(db_name)
    bucket_name = 'rewotes-pfu-bucket'
    
    @classmethod
    @abstractmethod
    def table_name(cls):
        """Name of the underlying SQL database table.
        :return: string: receiver's database table name
        """
        pass

    @classmethod
    @abstractmethod
    def create_table_sql(cls):
        """All the sql create scripts needed by file objects for tables and indices
        :return: string[] sql statements
        """
        pass
    
    @classmethod
    def create_tables(cls):
        """Creates receiver's table in the database if they don't already exist.
        :return: None
        """
        cursor = cls.db_conn.cursor()
        try:
            sqls = cls.create_table_sql()
            for sql in sqls: cursor.execute(sql)
            cls.db_conn.commit()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()

    @classmethod
    def run_sql_query(cls, sql):
        """Performs the given SQL query on the database and returns
           results as provided by the database.
        :param sql: SQL statement
        :return: tuple[]:  Array of tuples containing properties
        """
        cursor = cls.db_conn.cursor()
        try:
            return cursor.execute(sql).fetchall()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()

    @classmethod
    def run_sql_command(cls, sql):
        """Executes the given SQL on the database and commits.
        :param sql: SQL statement
        :return: None
        """
        cursor = cls.db_conn.cursor()
        try:
            cursor.execute(sql)
            cls.db_conn.commit()
        except sqlite3.Error as error:
            logger.error(f"Error running sql: {error}; ${sql}")
        finally:
            cursor.close()

    def __init__(self):
        pass
