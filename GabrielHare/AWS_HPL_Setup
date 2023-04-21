# AWS Parallel Cluster Setup
[Following general Setup instructions from AWS docs](https://docs.aws.amazon.com/parallelcluster/latest/ug/install-v3-virtual-environment.html)
Local Platform Only
Ubuntu 20.04 (Focal Fossa)

NOTE: Performance improvement option:  
[https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/efa.html#efa-instance-types](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/efa.html#efa-instance-types)  

NOTE: Logging option:  
[https://docs.aws.amazon.com/vpc/latest/userguide/flow-logs.html](https://docs.aws.amazon.com/vpc/latest/userguide/flow-logs.html)

NOTE: Relevant example use-cases:
[https://sunhwan.github.io/blog/2021/04/17/Run-Molecular-Dynamics-Simulation-on-AWS-Cluster.html](https://sunhwan.github.io/blog/2021/04/17/Run-Molecular-Dynamics-Simulation-on-AWS-Cluster.html)  
[https://aws.amazon.com/blogs/hpc/running-20k-simulations-in-3-days-with-aws-batch/](https://aws.amazon.com/blogs/hpc/running-20k-simulations-in-3-days-with-aws-batch/)  

## Install
(1) Install python3, pip, & virtualenv

    $ apt install python3
    $ python3 -m pip install --upgrade pip
    $ python3 -m pip install --user --upgrade virtualenv


(2) Install AWS Parallel Cluster into virtual environment apc

    $ mkdir [your cluster project] && cd [your cluster project]
    $ virtualenv apc
    $ source apc/bin/activate
    $ pip install --upgrade "aws-parallelcluster"
 
 Verify:

    $ pcluster version

  
(3) Install Node Version Manager  

    $ curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.38.0/install.sh | bash  
    $ chmod ug+x ~/.nvm/nvm.sh  
    $ source ~/.nvm/nvm.sh  
    $ nvm install --lts  

Verify:  

    $ node --version  

 
 ## Deploy cluster
(4) Configure & Create cluster  

    $ pcluster configure --config [your cluster config].yaml

This will start a command line interface to specify the cluster configuration.
EXAMPLE:
 - One queue, named "queue1"
 - Nodes at t2.micro instances
 - Platform is Ubuntu 20.04 (focal fossa)
 - Head node has public IP address, cluster nodes are private
 - Maximum of 16 instances

    $ pcluster create-cluster --cluster-configuration [your cluster config].yaml --cluster-name [your cluster name]

Expected response:
 
    {  
    "cluster": {  
    "clusterName": "[your cluster name]",  
    "cloudformationStackStatus": "CREATE_IN_PROGRESS",  
    "cloudformationStackArn": "arn:aws:cloudformation:us-west-1:918364550460:stack/mat3ra-benchmark/e06cb140-dbfa-11ed-8526-0658b83a89f5",  
    "region": "us-west-1",  
    "version": "3.5.1",  
    "clusterStatus": "CREATE_IN_PROGRESS",  
    "scheduler": {  
    "type": "slurm"  
    }  
    }  
    }  
Wait for cluster creation to complete

    $ watch -n 10 pcluster describe-cluster --cluster-name mat3ra-benchmark  

Cluster creation is complete when: 

    "clusterStatus": "CREATE_COMPLETE"  


NOTE: The EC2 console will show only the cluster head node. Additional nodes will be created as needed, and will be terminated after a period of inactivity.

NOTE: The configuration file will persist, and can be used to recreate the cluster, provided that the VPC identified in the configuration still exists.

## Log in
(5) Log in to cluster head node

    $ pcluster ssh --cluster-name mat3ra-benchmark

WARNING: DO NOT UPGRADE - this will stall out at 79%, and the cluster will need to be recreated!

# HPL Benchmark setup
[Following general setup instructions, with modification for AWS instances](https://alan-turing-institute.github.io/data-science-benchmarking/examples/HPL_benchmarks_workflow_example.html)
Remote Platform Only
Ubuntu 20.04 (focal fossa)

## Install linear algebra library
(1) Confirm [ATLAS](https://math-atlas.sourceforge.net) development library is installed

    $ dpkg --listfiles libatlas-base-dev | grep 'atlas_buildinfo\.h'

Expected response:

    /usr/include/x86_64-linux-gnu/atlas/atlas_buildinfo.h  

EXPECT: Install include path is:

    /usr/include/x86_64-linux-gnu/atlas/  

If ATLAS is not found install for development
$ apt install -y libatlas-base-dev  

## Install MPI library
(2) Confirm OpenMPI [LINK] is installed and active

    $ which mpirun

Expected response:

    /opt/amazon/openmpi/bin/mpirun

If [OpenMPI](https://docs.open-mpi.org) is not found [install for development](https://packages.ubuntu.com/focal/libopenmpi-dev)

    $ apt install -y openmpi-bin libopenmpi-dev

(2a) OPTIONAL: Enable Intel [MPI](https://www.intel.com/content/www/us/en/docs/mpi-library/get-started-guide-linux/2021-6/overview.html)
Confirm that IntelMPI is [available on cluster nodes](https://docs.aws.amazon.com/parallelcluster/latest/ug/intelmpi.html) as a [module](https://uisapp2.iu.edu/confluence-prd/pages/viewpage.action?pageId=115540061)

    $ module avail

Expect response to include "intelmpi"
Switch to using IntelMPI

    $ module load intelmpi

Confirm switch

    $ which mpirun

Expect response:

    /opt/intel/mpi/2021.6.0/bin/mpirun

IMPORTANT: IntelMPI will need to be enabled on each node as a part of the launch script, otherwise HPL will not be able to run.

NOTE: Another MPI option is [MPIch](https://www.mpich.org/about/overview/)
[https://stackoverflow.com/a/25493270/5917300](https://stackoverflow.com/a/25493270/5917300)

## Install HPL

[Following these directions for HPL install](https://gist.github.com/Levi-Hope/27b9c32cc5c9ded78fff3f155fc7b5ea)  
NOTE: "wget," "nano," and "build-essential" are pre-installed on AWS nodes

(3) Confirm latest HPL version in repostory (currently 2.3):
[https://netlib.org/benchmark/hpl/](https://netlib.org/benchmark/hpl/)  

(4) Download HPL source to user directory, which will be shared by all nodes.

    $ cd ~
    $ wget https://www.netlib.org/benchmark/hpl/hpl-2.3.tar.gz 
    $ tar -xf hpl-2.3.tar.gz

(4) Configure HPL build

    $ mv hpl-2.3 hpl  
    $ cd hpl/setup  
    $ sh make_generic  
    $ mv Make.UNKNOWN ../Make.linux  
    $ cd ../

Edit the Make.linux file
NOTE: In Make.linux "TOPdir = ${HOME}/hpl/", which is the reason for the rename above

    $ nano Make.linux

>     ARCH = linux  
>     LAinc = /usr/include/x86_64-linux-gnu/atlas/

If using Intel MPI  

>     MPinc = /opt/intel/mpi/2021.6.0/include/  
>     MPlib = /opt/intel/mpi/2021.6.0/lib/release/libmpi.so

If using OpenMPI  

>     MPinc = /opt/amazon/openmpu/include/  
>     MPlib = /opt/amazon/openmpu/lib/libmpi.so

Build xhpl  

    $ make arch=linux 

NOTE: "Expect gcc: warning: ..." but no errors
IMPORTANT: The path to "xhpl" is "${HOME}/hpl/bin/linux" so it must be added to path of each instance

Confirm build using a test for one node

    $ cd bin/linux
    $ mv HPL.dat default_HPL.dat
    $ nano HPL.dat

Paste the following into HPL.dat (every line must be included)

>     HPLinpack benchmark input file
>     Innovative Computing Laboratory, University of Tennessee
>     HPL.out      output file name (if any)
>     6            device out (6=stdout,7=stderr,file)
>     1            # of problems sizes (N)
>     5040         Ns
>     1            # of NBs
>     128          NBs
>     0            PMAP process mapping (0=Row-,1=Column-major)
>     1            # of process grids (P x Q)
>     1            Ps
>     1            Qs
>     16.0         threshold
>     1            # of panel fact
>     2            PFACTs (0=left, 1=Crout, 2=Right)
>     1            # of recursive stopping criterium
>     4            NBMINs (>= 1)
>     1            # of panels in recursion
>     2            NDIVs
>     1            # of recursive panel fact.
>     1            RFACTs (0=left, 1=Crout, 2=Right)
>     1            # of broadcast
>     1            BCASTs (0=1rg,1=1rM,2=2rg,3=2rM,4=Lng,5=LnM)
>     1            # of lookahead depth
>     1            DEPTHs (>=0)
>     2            SWAP (0=bin-exch,1=long,2=mix)
>     64           swapping threshold
>     0            L1 in (0=transposed,1=no-transposed) form
>     0            U  in (0=transposed,1=no-transposed) form
>     1            Equilibration (0=no,1=yes)
>     8            memory alignment in double (> 0)

Execute the test

    $ xhpl

  
NOTE: [Explanation of HPL.dat parameters](https://ulhpc-tutorials.readthedocs.io/en/latest/parallel/mpi/HPL/#hpl-main-parameters)
This includes instructions for choosing the maximum viable matrix size based on available memory in node hardware.

NOTE: [Configuring HPL for general node hardware](https://github.com/matthew-li/lbnl_hpl_doc#3-gathering-parameters)  

# Install Exabyte Benchmark Suite

NOTE: SLURM is already installed on AWS instances, and [no other systems resource management services are supported](https://aws.amazon.com/blogs/hpc/choosing-between-batch-or-parallelcluster-for-hpc/).
NOTE: AWS offers [a migration guide for older clusters](https://aws.amazon.com/blogs/hpc/easing-your-migration-from-sge-to-slurm-in-aws-parallelcluster-3/).
NOTE: SLURM provides a "qsub" command that accepts canonical PBS arguments.
If sbatch is not found

(1) Install git-lfs
[Following GitHub directions](https://github.com/git-lfs/git-lfs/blob/main/INSTALLING.md)
NOTE: Instances have git installed, but git-lfs is required for this project.

    $ curl -s [https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh](https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh) | sudo bash  
    $ sudo apt install -y git-lfs

(2) Install python3-virtualenv
NOTE: Instances will have python3 installed, but use of virtualenv is recommended for requirement installs

    $ sudo apt install -y python3-virtualenv

(3) Install project dependencies

    $ apt install -y libpng-dev libfreetype-dev  

NOTE: These libraries are already installed on instance

(4) Clone the benchmark repository
NOTE: This environment is not expected to be used for development, so a minimal clone will suffice.
TEMP: Use project from GabrielHare until update PR is merged

    $ git clone --depth 1 --recurse-submodules --shallow-submodules [https://github.com/GabrielHare/exabyte-benchmarks-suite.git](https://github.com/GabrielHare/exabyte-benchmarks-suite.git) --single-branch --branch GabrielHare
    $ cd exabyte-benchmarks-suite
    $ virtualenv env
    $ source env/bin/activate
    $ pip install -r requirements.txt

(5) Run the benchmark tests
NOTE: Only the HPL benchmarks will be run
 - espresso module is empty
 - vasp module expects additional files which are missing 
 - gromacs has not been installed

Prepare the benchmarks (generates HPL.dat files)

    $ ./exabench --prepare --type hpl  

Execute the tests (will launch EC2 node instances)

    $ ./exabench --execute --type hpl

Wait for queue to empty

    $ watch -n 10 squeue

Collect results when the queue is empty so all benchmark test configurations have executed.

    $ ./exabench --results --type hpl  

Results are appended to the results/result.json file, which is 

    $ git config diff.lfs.textconv cat  
    $ git diff results/results.json

In the default configuration EC2 node instances will terminate after 10 minutes without activity.

NOTE: Node instances will appear in EC2 to handle each job, up to the cluster maximum.
NOTE: Standard and error outputs will appear in a log file adjacent to the HPL.dat for each benchmark
NOTE: [SLURM squeue job status](https://slurm.schedmd.com/squeue.html#lbAG)
NOTE: [Additional diagnostics using sacctmgr](https://stackoverflow.com/questions/29928925/how-can-i-get-detailed-job-run-info-from-slurm-e-g-like-that-produced-for-sta)
NOTE: [Additional node boot failure information](https://stackoverflow.com/questions/59074208/obtain-the-boot-and-failure-history-of-nodes-in-a-slurm-cluster)
NOTE: The following is a *rare* job failure, and can be resolved simply be retrying the failed job:
```
--------------------------------------------------------------------------  
There are not enough slots available in the system to satisfy the 16  
slots that were requested by the application:  
  
xhpl  
  
Either request fewer slots for your application, or make more slots  
available for use.  
  
A "slot" is the Open MPI term for an allocatable unit where we can  
launch a process. The number of slots available are defined by the  
environment in which Open MPI processes are run:  
  
1. Hostfile, via "slots=N" clauses (N defaults to number of  
processor cores if not provided)  
2. The --host command line parameter, via a ":N" suffix on the  
hostname (N defaults to 1 if not provided)  
3. Resource manager (e.g., SLURM, PBS/Torque, LSF, etc.)  
4. If none of a hostfile, the --host command line parameter, or an  
RM is present, Open MPI defaults to the number of processor cores  
  
In all the above cases, if you want Open MPI to default to the number  
of hardware threads instead of the number of processor cores, use the  
--use-hwthread-cpus option.  
  
Alternatively, you can use the --oversubscribe option to ignore the  
number of available slots when deciding the number of processes to  
launch.  
--------------------------------------------------------------------------
```
