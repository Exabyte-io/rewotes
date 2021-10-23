#!/bin/bash

pytest test/ --cov=basistron

regimes="experimental calculated"

for regime in $regimes; do
    echo "${regime}"
    python -m basistron.app \
        --xyz_path "h2.xyz" \
        --regime "${regime}" \
        --property vibrational_frequency \
        --tolerance 0.5
    python -m basistron.app \
        --xyz_path "h2.xyz" \
        --regime "${regime}" \
        --property homo_lumo_gap \
        --tolerance 0.5
    python -m basistron.app \
        --xyz_path "ch4.xyz" \
        --regime "${regime}" \
        --property homo_lumo_gap \
        --tolerance 0.5
#    python -m basistron.app \
#        --xyz_path "ch4.xyz" \
#        --regime "${regime}" \
#        --property vibrational_frequency \
#        --value 4407 \
#        --tolerance 0.5
done
