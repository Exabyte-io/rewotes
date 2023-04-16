#!/usr/bin/env bash
echo 'Tear down infrastructure'
ansible-playbook ../ansible/playbooks/terraform-destroy.yml
