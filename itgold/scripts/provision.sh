#!/usr/bin/env bash
echo 'Provisioning environment: dev'
ansible-playbook -i ../ansible/inventory/dev.yml ../ansible/playbooks/terraform-apply.yml $1 $2 $3 $4

ansible-playbook -i ../ansible/inventory/dev.yml ../ansible/playbooks/terraform-destroy.yml $1 $2 $3 $4
