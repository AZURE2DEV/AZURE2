#!/bin/bash

# AZURE2 Docker Build Script
docker build --ulimit nofile=1024 -t azure2 -f ../packaging/docker/Dockerfile.azure2
