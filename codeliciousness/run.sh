#!/bin/bash

python -m basistron.app \
    --xyz_path h2.xyz \
    --target_property homo_lumo_gap \
    --reference_value 100.0
