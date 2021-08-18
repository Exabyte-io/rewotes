import os
from dataclasses import dataclass
from typing import List


@dataclass
class FileDesc:
    path: str
    filename: str


def list_files(directory_path: str) -> List[FileDesc]:
    all_files = []
    top_files = os.listdir(directory_path)

    for file in top_files:
        path = os.path.join(directory_path, file)
        if os.path.isdir(path):
            all_files.extend(list_files(path))
        else:
            all_files.append(FileDesc(path=path, filename=file))

    return all_files
