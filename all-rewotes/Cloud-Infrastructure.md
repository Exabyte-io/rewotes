# Cloud Infrastructure (HPC/DevOps)

> Ideal candidate: skilled HPC engineer versed in cloud, HPC, and DevOps

# Overview

The aim of this task is to create a CI/CD pipeline (github workflow) that includes (i) deploying cloud infrastructure for cluster compute, (ii) configuring it for running HPC application(s), and (iii) running benchmarks for a set of distributed memory calculations.

# Requirements

1. A working CI/CD pipeline - e.g. GitHub action - able to deploy and configure an HPC cluster
2. An automated workflow (using a configurable Github action) to benchmark one or more HPC application on one or more cloud instance type

# Expectations

- The application may be relatively simple - e.g. Linpack, this is focused more on infrastructure
- Clean workflow logic

# Timeline

We leave exact timing to the candidate. Should fit Within 5 days total.

# User story

As a user of this CI/CD pipeline I can:

- initiate tests for a specific number of scenarios: e.g. 2 nodes, 16 core per node
- select the instance type to be used 

# Notes

- Commit early and often

# Suggestions

We suggest:

- using AWS as the cloud provider
- using Exabench as the source of benchmarks: https://github.com/Exabyte-io/exabyte-benchmarks-suite
- using CentOS or similar as operating system
- using Terraform for infrastructure management
