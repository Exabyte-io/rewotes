# =============== Variables =============== #
# aws api access key
variable "aws_access_key" { type = string }

# aws api access secret
variable "aws_secret_key" { type = string }

# template to be used for ec2 instances
variable "instance_template_index" {
  type = number
  default = 1
}

# number of ec2 instances to instantiate
variable "instance_count" {
  type = number
  default = 1
}

terraform {
  required_version = ">= 1.0"
}

# initialize AWS provider
provider "aws" {
  region = "us-east-2"

  access_key = var.aws_access_key
  secret_key = var.aws_secret_key
}

# initialize compute networking and infrastructure
module "infra" {
  source = "./infra"
}

# initialize EC2 compute instances
module "compute" {
  source = "./compute"

  subnet_id = module.infra.ir_subnet_id
  security_group_id = module.infra.ir_security_group_id
  instance_template_index = var.instance_template_index
  instance_count = var.instance_count
}

# define output as an IP list of initialized instances
output "instance_ips" {
  value = module.compute.instance_ips
}
