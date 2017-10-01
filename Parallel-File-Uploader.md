# Multi-threaded file uploader (Backend)

> Ideal candidate: skilled python developer with solid knowledge of cloud and distributed systems.

# Overview

Create a python application that uploads a set of given files to a cloud object storage in parallel through the cloud provider's or third party API.

# Requirements

1. Support up to 100,000nds of files, all inside one directory with arbitrary sizes. The root directory may contain subdirectories.
1. The object storage container which holds the objects is private and only credential-based access is allowed.
1. Each object inside object storage should have an associated metadata which contains file size, last modification time and file permissions.

# Expectations

- Fast (utilize full network bandwidth), low CPU (do not block all other processes) and low Memory (<25% tentatively) file uploader
- Support for AWS S3
- Modular and Object oriented implementation (to add other cloud providers)
- Clean and documented code
- Tests

# Timeline

3-5 days total.

# Notes

We can provide temporary credentials to access AWS/Azure.
