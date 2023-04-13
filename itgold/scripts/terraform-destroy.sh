#!/usr/bin/env bash
echo 'Destroying infrastructure: dev'
export TF_VAR_aws_access_key
export TF_VAR_aws_secret_key
source .env
cd ../terraform || exit
terraform destroy -target=module.compute