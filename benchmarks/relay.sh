docker run --rm \
    -v $(pwd):/app \
    -w /app \
    omnipath-image:latest \
    python nest_relay.py
