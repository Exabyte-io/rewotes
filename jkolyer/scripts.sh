
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
    if [ $# -eq 0 ]; then
        echo "Enter directory tree depth:  ex. samplegen 3"
        exit 1
    fi
    flush_db
    rm -r samples && mkdir samples
    python jkolyer/directories.py --tree_depth $1 samples    
}

localstack_p() {
    if [ $# -eq 0 ]; then
        echo "Enter endpoint port:  ex. localstack_p 4566"
    else
        # run uploads with multiprocessing against localstack
        nice python pfu.py --root_dir ./samples --endpoint_url "http://localhost:$1" 
    fi
}

localstack_c() {
    if [ $# -eq 0 ]; then
        echo "Enter endpoint port:  ex. localstack_c 4566"
    else
        # run uploads with concurrency against localstack
        nice python pfu.py --root_dir ./samples --endpoint_url "http://localhost:$1" --concurrent
    fi
}

