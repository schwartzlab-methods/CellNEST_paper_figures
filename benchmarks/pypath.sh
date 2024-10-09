docker run --rm \
    -v $(pwd):/app \
    -w /app \
    pypath-image:latest \
    python pypath_relay.py
