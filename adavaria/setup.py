from setuptools import setup, find_packages
import os
from pathlib import Path


def package_files(directory, pattern='**/*'):
    paths = []
    for filepath in Path(directory).glob(pattern):
        if filepath.is_file():
            paths.append(str(filepath.relative_to(directory)))
    return paths

setup(
    name='mlband',
    version=0.01,
    packages=find_packages(),
    package_data={
        'mlband': ['atom_init.json', 'magpie_data.csv'],
        'mlband': ['pretrained_model/*'],
    },
    author='Ali Davariashtiyani',
    author_email='a.davari871@gmail.com',
    description='CGCNN for band gap prediction',
    python_requires='>=3.7',
)