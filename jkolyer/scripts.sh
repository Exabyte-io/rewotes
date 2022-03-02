
runtest() {
    pytest -s --log-cli-level=DEBUG --asyncio-mode=auto test
}

samplegen() {
    rm -r samples && mkdir samples
    python jkolyer/directories.py --tree_depth 3 samples    
}

locakstack() {
    python pfu.py --root_dir ./samples --endpoint_url "http://localhost:4566"
}

