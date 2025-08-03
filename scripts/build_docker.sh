#!/bin/bash

# AZURE2 Docker Build Script
docker build --ulimit nofile=1024 -t azure2 -f ../docker/Dockerfile.azure2
