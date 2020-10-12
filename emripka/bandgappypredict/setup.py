import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    #name="bandgappypredict-YOUR-USERNAME-HERE", # Replace with your own username
    name="bandgappypredict", 
    version="0.0.1",
    author="EmilyRipka",
    author_email="emripka@gmail.com",
    description="A package to predict a material's bandgap for a given set of stoichiometric and physical parameters.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
