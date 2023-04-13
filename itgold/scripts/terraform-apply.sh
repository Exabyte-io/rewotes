#!/usr/bin/env bash
echo 'Applying terraform: dev'
export TF_VAR_aws_access_key
export TF_VAR_aws_secret_key
source .env
cd ../terraform || exit
terraform apply