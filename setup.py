# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

# The directory containing this file
HERE = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# This call to setup() does all the work
setup(
    name="mlbands",
    version="1.0.0",
    description="A Python package that implements automatic prediction of electronic band gaps for a set of materials based on training data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/andrewrgarcia/rewotes/tree/andrewrgarcia",
    author="Andrew Garcia, PhD",
    license="MIT",
    classifiers=[
        "Intended Audience :: Information Technology",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Operating System :: OS Independent"
    ],
    packages=["mlbands"],
    include_package_data=True,
    install_requires=["numpy"]
)
