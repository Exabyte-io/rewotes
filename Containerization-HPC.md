# Containerization / Benchmarks (HPC)

> Ideal candidate: skilled HPC engineer versed in HPC, and Containers

# Overview

The aim of this task is to build an HPC compatible container (i.e. [Singularity](https://sylabs.io/guides/3.5/user-guide/introduction.html)) and test its performance in comparison with a native installation (no containerization) for a set of distributed memory calculations.

# Solution

Here we use AWS virtual parallel cluster (https://docs.aws.amazon.com/parallelcluster/latest/ug/what-is-aws-parallelcluster.html). AWS account is needed. We suppose to use Linux and console tools. aws console tools are supposed to be installed.

**Note**: aws account should have appropriate permissions (TO BE ADDED).

## Steps

1. Install pcluster tool (here v2 is used): `pip3 install "aws-parallelcluster<3.0" --upgrade --user`
2. Run initial conigure: `aws configure; pcluster configure`. The scheduler preferred to be slurm
3. Create a cluster: `pcluster create my-hpc-cluster`
4. Go to the cluster `pcluster ssh my-hpc-cluster`
5. Install the singularity from github releases (https://github.com/sylabs/singularity/releases/): `wget https://github.com/sylabs/singularity/releases/download/v3.9.8/singularity-ce_3.9.8-focal_amd64.deb; sudo apt install singularity-ce_3.9.8-focal_amd64.deb`
6. Copy and edit reciepe file (look the [hpl1.def](hpl1.def))
7. Build a singularity image: `sudo singularity build my-image.sif my-reciepe.def` (note **sudo** usage).
8. Run the application like this: `sbatch -n NUM_CPUS ./sing-run.sh --bind INPUT:/INPUT my-image.sif /path/to/app/in/image`. Here `--bind ...` is optional, it passes the input file into the container, can be specified multiple times and can pass directories.

`sing-run.sh` can be modified for custom needs.

**Important note**: do not use standard openmpi, you should use amazon toolkit and amazon openmpi. See the example reciepe file (hpl1.def) to details.

## Options

You can prepare your singularity image on another aws instance or computer, but you should optimize it for aws instance processor (AVX512 preferred).

You can set up the cluster without slurm and run the app via mpirun directly.

You can set up the cluster without slurm via terraform and salt/ansible, but details should be specified in later investigations.

## Example for HPL

Here is commented hpl1.def file, use it for new images as a sample.

```
Bootstrap: docker
# use amazon linux, it is based on centos. If use ubuntu or alpine, change packages management.
From: amazonlinux

# files to be copied into container
%files
  hpl-2.3.tar.gz
  Make.t1
 
# default environment, add amazon openmpi path
%environment
  export LC_ALL=C
  export PATH=$PATH:/opt/amazon/openmpi/bin

# prepare our image
%post
  export PATH=$PATH:/opt/amazon/openmpi/bin
  
  # install gcc etc
  yum groupinstall 'Development Tools' -y
  
  # install additional tools and atlas lib (for HPL)
  yum install -y make tar curl sudo altas atlas-devel
  
  # Here is ubuntu option
  # apt install -y altas make tar openmpi-aws
  
  # get and install amazon efa tools. Important for MPI applications
  curl -O https://efa-installer.amazonaws.com/aws-efa-installer-1.15.1.tar.gz
  tar -xf aws-efa-installer-1.15.1.tar.gz && pushd aws-efa-installer
  sudo ./efa_installer.sh -y
  popd
  
  #
  # Here is HPL-oriented part. Ypu can use another commands to build and install your app
  # build HPL in the /hpl directory, then we'll delete it
  mkdir /hpl
  tar xfz hpl-2.3.tar.gz -C /hpl
  rm hpl-2.3.tar.gz
  cp Make.t1 /hpl/hpl-2.3
  cd /hpl/hpl-2.3
  
  # use cutom prepared makefile
  make arch=t1
  
  # copy compiled xhpl app
  mv bin/t1/xhpl /bin
  cd ..
  
  # clean
  rm -rf /hpl
```

## Testing results

See [singularity-hpl.out](singularity-hpl.out) and [raw-hpl.out](raw-hpl.out) for singularity run and raw hpl run respectively. The short data is here:

```
RAW HPL:

WR10L2L2        5000   512     9     8               1.73             4.8320e+01
WR10L2L2       20000   512     9     8              22.08             2.4160e+02
WR10L2L2       12128   512     9     8               7.07             1.6818e+02

Singularity HPL:

WR10L2L2        5000   512     9     8               1.79             4.6495e+01
WR10L2L2       20000   512     9     8              22.52             2.3682e+02
WR10L2L2       12128   512     9     8               7.09             1.6776e+02
```

## Terraform approach notes

Here is my attempts to set up the cluster. See my terraform config [...](...). Please, note, that amazon efa tools should me installen on each node! Salt config should be added.

## Attached files

- [hpl1.def](hpl1.def) - sigularity reciepe to build an image
- [HPL.dat](HPL.dat) - Linpack input sample data
- [Make.t1](Make.t1) - makefile for Linpack (note path to atlas, in case of ubuntu it will be different)
- [raw-hpl.out](raw-hpl.out) - output of native HPL run
- [singularity-hpl.out](singularity-hpl.out) - output of singularity HPL run
- [results.txt](results.txt) - short results for comparation (native and singularity HPL)
- [xhpl.sh](xhpl.sh) - batch script for SLURM to run singilarity xhpl
- [zhpl.sh](zhpl.sh) - batch script for SLURM to run native xhpl
- [sing-run.sh](sing-run.sh) - batch script for SLURM to run any singularity container (just add singularity options, needed after `exec`).

