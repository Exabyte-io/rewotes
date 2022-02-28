
runtest() {
    pytest -s  --log-cli-level=DEBUG test
}

samplegen() {
    rm -r samples && mkdir samples
    python jkolyer/directories.py --tree_depth 3 samples    
}
