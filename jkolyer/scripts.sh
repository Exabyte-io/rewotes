
activate() {
    source bin/activate    
}

runtest() {
    pytest -s --log-cli-level=DEBUG --asyncio-mode=auto test
}

flush_db() {
    rm -r parallel-file-upload.db
}

samplegen() {
    flush_db
    rm -r samples && mkdir samples
    python jkolyer/directories.py --tree_depth $1 samples    
}

localstack_p() {
    # run uploads with multiprocessing against localstack
    python pfu.py --root_dir ./samples --endpoint_url "http://localhost:$1" 
}

localstack_c() {
    # run uploads with concurrency against localstack
    python pfu.py --root_dir ./samples --endpoint_url "http://localhost:$1" --concurrent
}

