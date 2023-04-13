# =============== Variables =============== #
# aws api access key
variable "aws_access_key" { type = string }

# aws api access secret
variable "aws_secret_key" { type = string }

terraform {
  required_version = ">= 1.0"
}

# initialize AWS provider
provider "aws" {
  region = "us-east-2"

  access_key = var.aws_access_key
  secret_key = var.aws_secret_key
}


module "infra" {
  source = "./infra"
}

module "compute" {
  source = "./compute"

  subnet_id = module.infra.ir_subnet_id
  instance_template_number = 1
}