BasisTron
=========

TODO
----

* Basis
  - how to update workflow model with basis set info

* App
  - default reference level of theory, basis set
  - compute allowed basis sets from reference data
  - return job id of submitted job
  - extras
    - match ranked basis sets to allowed basis sets
    - pick most compact allowed set from ranked set

BasisTron is the automatic basis set selection tool you've
always needed but have never had the time to write yourself.

Since this is a POC rather than a robust operational program,
please follow these manual setup instructions to ensure your
BasisTron experience is a smooth one. BasisTron relies on two
databases, a basis set database and a reference data
database. For the purposes of this project, the basis set
database is provided by the EMSL Basis Set Exchange and the
reference data "database" is provided by CCCBDB.

Installing the basis set database
---------------------------------

* Navigate to https://www.basissetexchange.org
* Click on the Download button at the top of the page
  - Choose NWChem as the Basis Set Format
  - Choose tar + bz2 as the Archive Type
* Press the Download button

After you have downloaded the basis set tarball, follow these
steps in a terminal (assuming your tarball was downloaded to
`~/Downloads`).

```bash
VERSION=v0.8.13 # at the time of project inception
mkdir -p ~/.basistron/basis/
mv ~/Downloads/basis_sets-nwchem-${VERSION}.tar.bz2 ~/.basistron/basis/
cd ~/.basistron/basis/
bunzip2 basis_sets-nwchem-${VERSION}.tar.bz2
tar -xvf basis_sets-nwchem-${VERSION}.tar
```

This provides the basis set database that BasisTron uses to
systematically rank and choose basis sets for a given system.


Installing the reference data database
--------------------------------------

There is no conveniently obtained "dump" for the CCCBDB. Therefore,
a small API client serves to dynamically fetch relevant reference
data at run-time. Query results are persisted to disk so subsequent
calls for the same data avoid the network. Persisted queries are
stored in `~/.basistron/cccbdb/` internally and should not be
accessed outside of the provided API.


Program Usage
=============

The environment for this program is managed with `poetry`. It can
be installed using pip into a matching python version.

```bash
$ python --version # ensure python in your path is ~3.9
Python 3.9.x
$ python -m pip install poetry
...
$ poetry install
...
$ poetry shell
```

The `poetry shell` command spawns a subshell with a virtualenv-like
experience. Then the BasisTron program can be executed using the
following command pattern:

```bash
export EXABYTE_USERNAME=yourusername
export EXABYTE_PASSWORD=yourpassword

python -m basistron.app \
    --xyz_path /path/to/file \
    --target_property homo_lumo_gap \
    --reference_value 2.0
```

The program assumes the contents of the XYZ file are in units of
angstroms and constitute a neutral singlet electronic configuration.

