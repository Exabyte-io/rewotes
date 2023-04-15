#!/usr/bin/env bash
echo 'Checkout projects: local'
ansible-playbook -i ../ansible/inventory/dev.yml ../ansible/playbooks/terraform-apply.yml $1 $2 $3 $4
ansible-playbook -i ../ansible/inventory/dev.yml ../ansible/playbooks/provision.yml $1 $2 $3 $4
ansible-playbook -i ../ansible/inventory/dev.yml ../ansible/playbooks/quantum_espresso.yml $1 $2 $3 $4
#ansible-playbook -i ../ansible/inventory/dev.yml ../ansible/playbooks/exabyte.yml $1 $2 $3 $4
#ansible-playbook -i ../ansible/inventory/dev.yml ../ansible/playbooks/terraform-destroy.yml $1 $2 $3 $4
