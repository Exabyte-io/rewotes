# Instructions

BasisSetSelector can be used as python module, or as command line script `optimize_basis`. To install, run 

```
pip install . 
```
in the root directory (i.e. the same directory as pyproject.toml).

Setting up the calculation can be done using a properly formatted .json. These examples assume NWChem is installed on your machine.

Working examples:

```bash
  # python script
  python example.py
  # command line script
  optimize_basis ex1.json b3lyp 0.1 --verbose
```
