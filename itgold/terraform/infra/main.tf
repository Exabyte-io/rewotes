# =============== Variables =============== #

# =============== AWS infrastructure =============== #
# Create IR VPC
resource "aws_vpc" "ir_vpc" {
  cidr_block = "10.0.0.0/16"

  tags = {
    Name = "ir_vpc"
  }
}

# Create a security group
resource "aws_security_group" "ir_security_group" {
  name   = "ir_security_group"
  vpc_id = aws_vpc.ir_vpc.id

  egress {
    from_port   = 0
    to_port     = 65535
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }

  ingress {
    from_port   = 0
    to_port     = 65535
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = {
    Name = "ir_security_group"
  }
}

# Create network ACL
resource "aws_network_acl" "ir_network_acl" {
  vpc_id = aws_vpc.ir_vpc.id

  ingress {
    rule_no    = 100
    protocol   = "tcp"
    action     = "allow"
    cidr_block = "0.0.0.0/0"
    from_port  = 0
    to_port    = 65535
  }

  egress {
    rule_no    = 200
    protocol   = "all"
    action     = "allow"
    cidr_block = "0.0.0.0/0"
    from_port  = 0
    to_port    = 0
  }


  tags = {
    Name = "ir_network_acl"
  }
}

# Create a private subnet
resource "aws_subnet" "ir_subnet" {
  vpc_id     = aws_vpc.ir_vpc.id
  cidr_block = "10.0.1.0/24"

  tags = {
    Name = "ir_subnet"
  }
}

# Associate subnet with the network ACL
resource "aws_network_acl_association" "ir_network_acl_subnet_association" {
  subnet_id      = aws_subnet.ir_subnet.id
  network_acl_id = aws_network_acl.ir_network_acl.id
}

# Create network interface
resource "aws_network_interface" "ir_network_interface" {
  subnet_id       = aws_subnet.ir_subnet.id
  security_groups = [aws_security_group.ir_security_group.id]

  tags = {
    Name = "ir_network_interface"
  }
}

# Create VPC internet gateway
resource "aws_internet_gateway" "ir_internet_gateway" {
  vpc_id = aws_vpc.ir_vpc.id

  tags = {
    Name = "ir_internet_gateway"
  }
}

# Create route table for internet gateway
resource "aws_route_table" "ir_route_table" {
  vpc_id = aws_vpc.ir_vpc.id

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.ir_internet_gateway.id
  }

  tags = {
    Name = "ir_route_table"
  }
}

# Associate route table with the subnet
resource "aws_route_table_association" "ir_route_table_subnet_association" {
  subnet_id      = aws_subnet.ir_subnet.id
  route_table_id = aws_route_table.ir_route_table.id
}

output "ir_subnet_id" {
  description = "ir_subnet object id output"
  value       = aws_subnet.ir_subnet.id
}

output "ir_security_group_id" {
  description = "ir_security_group object id output"
  value       = aws_security_group.ir_security_group.id
}
