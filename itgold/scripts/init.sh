#!/usr/bin/env bash
echo "Usage: $0 <instance_template_index> <instance_count>" >&2

# Default values
instance_template_index=0
instance_count=1

# Check if parameters are passed
if [ $# -ge 1 ]; then
  instance_template_index=$1
else
  instance_template_index=1
fi

if [ $# -ge 2 ]; then
  instance_count=$2
else
  instance_count=1
fi

# Confirm the values
echo "Instance Template Index: ${instance_template_index}"
echo "Instance Count: ${instance_count}"

echo 'Starting AWS deployment pipeline'
ansible-playbook ../ansible/playbooks/terraform-create.yml --extra-vars "instance_template_index=${instance_template_index} instance_count=${instance_count}"

echo 'Collecting hosts to provision'
ansible-playbook ../ansible/playbooks/terraform-collect.yml

echo 'Provisioning compute hosts'
ansible-playbook ../ansible/playbooks/provision.yml

echo 'Installing Quantum Espresso'
ansible-playbook ../ansible/playbooks/quantum_espresso.yml

echo 'Installing Exabyte Benchmark Suite'
ansible-playbook ../ansible/playbooks/exabyte.yml

#echo 'Tear down infrastructure'
#ansible-playbook ../ansible/playbooks/terraform-destroy.yml
