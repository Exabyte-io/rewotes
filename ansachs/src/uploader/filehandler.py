import os
from dataclasses import dataclass


@dataclass
class FileDesc:
    path: str
    filename: str


def list_files(directory_path: str) -> list[FileDesc]:
    all_files = list()
    top_files = os.listdir(directory_path)

    for file in top_files:
        path = os.path.join(directory_path, file)
        if os.path.isdir(path):
            all_files.extend(list_files(path))
        else:
            all_files.append(FileDesc(path=path, filename=file))

    return all_files
