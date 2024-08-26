FROM python:3.9-slim
WORKDIR /app

RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

# forked version to accomodate NEST's custom db 
RUN pip install -U git+https://github.com/d-paliwal/SpatialDM.git 
RUN pip install pyyaml

COPY . /app
