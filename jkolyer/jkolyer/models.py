import cuid

class BaseModel:
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
                    fileSize INTEGER,
                    lastModified INTEGER,
                    permissions TEXT,
                    path TEXT,
                    fileName TEXT,
                    filePath TEXT
                  );
                """.format(table_name=cls.table_name()),
                'CREATE INDEX IF NOT EXISTS IdxFileName ON FileStat(fileName)',
                'CREATE INDEX IF NOT EXISTS IdxFilePath ON FileStat(filePath)',
                'CREATE INDEX IF NOT EXISTS IdxFileSize ON FileStat(fileSize)'
                ]

class JobModel(BaseModel):
    @classmethod
    def table_name(cls):
        return 'UploadJob'
    
    @classmethod
    def create_table_sql(cls):
        return ["""
                  CREATE TABLE IF NOT EXISTS {table_name}
                  ( id TEXT PRIMARY KEY, 
                    fileId TEXT,
                    status TEXT,
                    createdAt INTEGER,
                    updatedAt INTEGER,
                    FOREIGN KEY (fileId) REFERENCES {file_table_name}(id)
                  );
               """.format(table_name=cls.table_name(), file_table_name=FileModel.table_name()),
                'CREATE INDEX IF NOT EXISTS IdxJobFile ON UploadJob(fileId);',
                'CREATE INDEX IF NOT EXISTS IdxStatus ON UploadJob(status);'
                ]
