# =============== Variables =============== #
variable "subnet_id" {
}

variable "instance_template_number" {
  type    = number
}

# types of ec2 instances declaration
variable "instance_templates" {
  description = "A map of EC2 instances"
  type        = map(object({
    ami           = string
    instance_type = string
  }))

  # todo: move to configuration
  default = {
    "1" = {
      ami           = "ami-0a695f0d95cefc163"
      instance_type = "t2.micro"
    },

    "2" = {
      ami           = "ami-0103f211a154d64a6"
      instance_type = "t2.micro"
    },

    "3" = {
      ami           = "ami-0fb3a91b7ce257ec1"
      instance_type = "t2.micro"
    }
  }
}

# Defines a path to ec2 hosts access key
variable "key_path" {
  default = "access/ir-keypair.pem"
}

output "instance_template_number" {
  description = "instance template number output"
  value       = var.instance_template_number
}

# =============== EC2 instances =============== #

# Create a simple EC2 instance
resource "aws_instance" "ir_ec2_instance" {
  ami           = var.instance_templates[var.instance_template_number].ami
  instance_type = var.instance_templates[var.instance_template_number].instance_type
  subnet_id     = var.subnet_id
  count         = 1

  connection {
    type        = "ssh"
    user        = "ubuntu"
    private_key = file(var.key_path)
    host        = self.public_ip
  }

  tags = {
    Name = "ir_ec2_instance-${count.index + 1}"
  }
}