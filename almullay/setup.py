# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

# with open('README.rst') as f:
#     readme = f.read()

# with open('LICENSE') as f:
#     license = f.read()

setup(
    name='sample',
    # version='0.1.0',
    # description='Sample',
    # long_description=readme,
    # author='Yousif Almulla',
    # author_email='yousif.almulla38@gmail.com',
    # url='',
    # license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)