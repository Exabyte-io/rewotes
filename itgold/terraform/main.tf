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
  security_group_id = module.infra.ir_security_group_id
  instance_template_number = 1
}

# populate list of public IPs for EC2 resources
resource "null_resource" "collect_ips" {
  depends_on = [module.compute.instance_ips]

  provisioner "local-exec" {
    command = <<EOF
      #!/bin/bash
      sudo apt-get update
      sudo apt-get install -y jq
      rm -rf ../compute.list || true
      terraform output -json instance_ips | jq -r '.[]' > ../compute.list
    EOF
  }
}
