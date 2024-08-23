FROM python:3.9-slim
WORKDIR /app

RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

RUN pip install -U git+https://github.com/StatBiomed/SpatialDM

COPY . /app
