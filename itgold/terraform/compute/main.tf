# =============== Variables =============== #
variable "subnet_id" {
}

variable "security_group_id" {
}

variable "instance_template_index" {
  type = number
  default = 1
}

variable "instance_count" {
  type = number
  default = 1
}

# types of ec2 instances declaration
variable "instance_templates" {
  description = "A map of EC2 instances"
  type        = map(object({
    ami            = string
    instance_type  = string
    cpu_core_count = number
  }))

  # todo: move to configuration
  default = {
    "1" = {
      ami            = "ami-0a695f0d95cefc163"
      instance_type  = "t2.micro"
      cpu_core_count = 1
    },

    "2" = {
      ami            = "ami-0a695f0d95cefc163"
      instance_type = "t2.micro"
      cpu_core_count = 1
    },

    "3" = {
      ami            = "ami-0a695f0d95cefc163"
      instance_type = "t2.micro"
      cpu_core_count = 1
    }
  }
}

# Defines ssh access key pair
resource "aws_key_pair" "ir_ssh_key_pair" {
  key_name   = "ir-key-pair"
  public_key = file("../access/id_rsa.pub")
}

output "instance_template_index" {
  description = "instance template number output"
  value       = var.instance_template_index
}

# =============== EC2 instances =============== #

# Create a simple EC2 instance
resource "aws_instance" "ir_ec2_instance" {
  ami                         = var.instance_templates[var.instance_template_index].ami
  instance_type               = var.instance_templates[var.instance_template_index].instance_type
  subnet_id                   = var.subnet_id
  key_name                    = aws_key_pair.ir_ssh_key_pair.key_name
  associate_public_ip_address = true
  count                       = var.instance_count

  vpc_security_group_ids = [var.security_group_id]

  connection {
    type        = "ssh"
    user        = "ubuntu"
    private_key = file("../access/id_rsa")
    host        = self.public_ip
  }
  tags = {
    Name = "ir_ec2_instance-${count.index + 1}"
  }

  # Creating user luke with the right permissions and ssh access
  provisioner "remote-exec" {
    inline = [
      "sudo useradd -m -s /bin/bash luke",
      "echo \"luke ALL=(ALL) NOPASSWD: ALL\" | sudo tee -a /etc/sudoers",
      "sudo mkdir -p /home/luke/.ssh",
      "sudo sh -c 'echo \"${file("../access/id_rsa.pub")}\" >> /home/luke/.ssh/authorized_keys'",
      "sudo chmod 700 /home/luke/.ssh",
      "sudo chmod 600 /home/luke/.ssh/authorized_keys",
      "sudo chown -R luke:luke /home/luke/.ssh"
    ]
  }
}

output "instance_ips" {
  value = aws_instance.ir_ec2_instance.*.public_ip
}
