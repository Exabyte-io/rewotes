# Cloud Infrastructure (HPC/DevOps)

# Overview

The aim of this task is to create a CI/CD pipeline (github workflow) that includes (i) deploying cloud infrastructure for cluster compute, (ii) configuring it for running HPC application(s), and (iii) running benchmarks for a set of distributed memory calculations.

# Requirements

1. A working CI/CD pipeline - e.g. GitHub action - able to deploy and configure an HPC cluster
2. An automated workflow (using a configurable Github action) to benchmark one or more HPC application on one or more cloud instance type

# User story

As a user of this CI/CD pipeline I can:

- initiate tests for a specific number of scenarios: e.g. 2 nodes, 16 core per node
- select the instance type to be used

# The Flow
1. When user pushes anything to the repository, it triggers github pipeline
2. Pipeline initiates AWS cluster with terraform
3. Pipeline initialize ec2 hosts
4. Pipeline deploys application Quantum Espresso
4. Pipeline deploys application Exabyte Test Suite
7. Pipeline shuts down the whole AWS infrastructure

# Project assumptions and allowances
- AWS infrastructure is not the most secure configuration. Some work needs to be done to make it better on that front
- The pipeline just builds and Quantum Espresso on all declared instances
- Pipeline builds Exabyte Test Suite but is not running it due to some libraries compatibility errors
- EC2 instance templates should be moved to configuration going forward 

# Project structure

_scripts_ - contains all the required scripts to run the pipeline:

_ansible_ - ansble playbooks. All pipeline steps are managed by ansible

_terraform_ - terraform scripts. Terraform is used to initialize AWS infrastructure and EC2 compute instances

_access_ - private/public keys to access EC2 instances

# Initialization
Deployment pipeline can be executed from GitHub Actions workflow as well as from the local environment command line

To run application from github there are 2 **Repository secrets** required to be initialized:

_AWS_ACCESS_KEY_ID - aws access key_

_AWS_SECRET_ACCESS_KEY - aws secret key_

It could be done by the following URL: https://github.com/itgold/rewotes/settings/secrets/actions by clicking on New Repository Secret.

To run it from local, **.env** file needs to be created in the _/scripts_ folder in the following format:

_TF_VAR_aws_access_key: aws_access_key_

_TF_VAR_aws_secret_key: aws_secret_key_

# Execution
_init.sh_ - perform most of the steps from The Flow except shutting down the infrastructure

_destroy.sh_ - performs AWS infrastructure shutting down operation

# Run example
Run _scripts/init.sh_ with 2 parameters

_instance_template_index_: what type of EC2 instance to run. List of available types is defined in _itgold/terraform/compute/main.tf_ file under instance_templates variable.  

_instance_count_: how many instances to instantiate

Example: 

_./init.sh 1 2_

Will execute the pipeline on ami-0a695f0d95cefc163, t2.micro, 1 cpu core. The pipeline will run on 2 hosts

And when you are ready to delete the whole infrastructure:

_./destroy.sh_

