#!/bin/bash
# ============================================================
# AZURE2 Docker Build Script for macOS
# Automatically starts Docker Desktop if needed and waits until it's ready
# Usage: bash build_docker.sh
# ============================================================

set -e  # Exit immediately if a command fails

echo "=== AZURE2 Docker Build Script ==="

# --- Check Docker installation ---
if ! command -v docker >/dev/null 2>&1; then
    echo "Docker is not installed. Please install Docker Desktop from:\n"
    echo "https://www.docker.com/products/docker-desktop \n"
    exit 1
fi

# --- Start Docker Desktop if not running ---
if ! docker info >/dev/null 2>&1; then
    echo "Docker is not running. Starting Docker Desktop... \n"
    open --background -a Docker

    echo -n "Waiting for Docker to start \n"
    until docker info >/dev/null 2>&1; do
        echo -n "."
        sleep 2
    done
    echo " Docker is ready. \n"
else
    echo "Docker is already running. \n"
fi

# --- Build the Docker image ---
echo "🔨 Building Docker image: azure2 \n"
docker build --ulimit nofile=1024 -t azure2 -f ../docker/Dockerfile.azure2 ..

if [ $? -eq 0 ]; then
    echo "Docker image 'azure2' built successfully! \n"
else
    echo "Docker build failed. Check the logs above for details.\n"
    exit 1
fi
