# Containerization / Benchmarks (HPC)

> Ideal candidate: skilled HPC engineer versed in HPC, and Containers

# Overview

The aim of this task is to build an HPC compatible container (i.e. [Singularity](https://sylabs.io/guides/3.5/user-guide/introduction.html)) and test its performance in comparison with a native installation (no containerization) for a set of distributed memory calculations.

# Requirements

1. A working deployment pipeline - using any preferred tool such as SaltStack, Terraform, CloudFormation - for building out the computational infrastructure
2. A pipeline for building the HPC compatible container
3. A set of benchmarks for one or more HPC application on one or more cloud instance type

# Expectations

- The application may be relatively simple - e.g. Linpack, this is focused more on infrastructure
- Repeatable approach (no manual setup "in console")
- Clean workflow logic

# Timeline

We leave exact timing to the candidate. Should fit Within 5 days total.

# User story

As a user of this pipeline I can:

- build an HPC-compatible container for an HPC executable/code
- run test calculations to assert working state of this container
- (optional) compare the behavior of this container with a OS native installation

# Notes

- Commit early and often

# Suggestions

We suggest:

- using AWS as the cloud provider
- using Exabench as the source of benchmarks: https://github.com/Exabyte-io/exabyte-benchmarks-suite
- using CentOS or similar as operating system
- using Salstack, or Terraform, for infrastructure management
