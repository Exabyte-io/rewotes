#!/usr/bin/env bash
echo 'Provisioning environment: dev'
ansible-playbook -i ../ansible/inventory/dev.yml ../ansible/playbooks/provision_host.yml $1 $2 $3 $4
