#!/bin/bash

python -m basistron.app \
    --xyz_path h2.xyz \
    --property homo_lumo_gap \
    --tolerance 20 \
    --reference_theory 'CCSD(T)' \
    --target_theory B3LYP \
    --regime calculated
