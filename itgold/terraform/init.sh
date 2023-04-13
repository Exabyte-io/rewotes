export TF_VAR_aws_access_key
export TF_VAR_aws_secret_key
source .env
aws ec2 create-key-pair --key-name ir-keypair --query 'KeyMaterial' --output text > access/ir-keypair.pem