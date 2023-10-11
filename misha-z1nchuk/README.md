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

We leave exact timing to the candidate. Must fit Within 5 days total.

# Notes

- we can provide temporary credentials to access AWS/Azure.



# Multithreaded File Uploader (AWS S3)




## Environment Variables

To run this project, you will need to add the following environment variables to your .env file

`AWS_BUCKET_NAME` `AWS_ACCESS_KEY` `AWS_ACCESS_SECRET_KEY`   `AWS_BUCKET_REGION`

See .env.example file.



## Run Locally

Clone the project

```bash
  git clone https://github.com/misha-z1nchuk/multithread-s3-file-uploader
```

Go to the project directory

```bash
  cd multithread-s3-file-uploader
```

Install dependencies

```bash
  npm install
```

Add folder that you want to upload to root folder of project

And run

```bash
    node lib/main.js -f ../<folder name>
```
