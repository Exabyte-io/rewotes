import setuptools

with open("README.md", "r") as readme:
    long_description = readme.read()

setuptools.setup(name="convtrack",
                 version="0.1a1",
                 author="James Dean",
                 description="Basic implementation of convergence with respect to supercell size",
                 long_description = long_description,
                 long_description_content_type = "test/markdown",
                 platforms=["Unix"],
                 packages=setuptools.find_packages(exclude="test"),
                 include_package_data=False,
                 keywords="chemistry",
                 classifiers=["Development Status :: 2 - Alpha",
                              "Environment :: Console",
                              "Intended Audience :: Science/Research",
                              "Operating System :: Unix",
                              "Programming Language :: Python :: 3.6",
                              "Topic :: Scientific/Engineering :: Chemistry"],
                 python_requires=">=3.6",
                 install_requires=["ase", "numpy"])