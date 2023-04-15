#!/usr/bin/env bash
echo 'Checkout projects: local'
#ansible-playbook ../ansible/playbooks/terraform-apply.yml $1 $2 $3 $4
ansible-playbook ../ansible/playbooks/terraform-collect.yml $1 $2 $3 $4
#ansible-playbook ../ansible/playbooks/provision.yml $1 $2 $3 $4
#ansible-playbook ../ansible/playbooks/quantum_espresso.yml $1 $2 $3 $4
#ansible-playbook ../ansible/playbooks/exabyte.yml $1 $2 $3 $4
#ansible-playbook ../ansible/playbooks/terraform-destroy.yml $1 $2 $3 $4
