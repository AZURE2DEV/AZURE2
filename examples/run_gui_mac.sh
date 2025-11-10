#!/bin/bash
# -------------------------------------------------------------
# Script to run AZURE2 GUI application inside Docker on macOS
# with XQuartz display forwarding.
#
# Usage:
#     ./run_gui_mac.sh
#
# Prerequisites:
#   • XQuartz installed and "Allow connections from network clients" enabled.
#   • Docker Desktop installed and running.
#   • GIT repository cloned:
#       git clone https://github.com/akhilb5/AZURE2.git
#       git checkout api
#    
#   • AZURE2 Docker image built:
#       docker build -t azure2:latest ../docker
# -------------------------------------------------------------

set -e

# --- Configurable Variables ---
CONTAINER_NAME=azure2_container # name of the docker container
HOST_OUTPUT=~/Desktop/azure2_output   # where outputs appear on host computer the mac 
APP_DIR="$PWD/azure2"                 # where your source lives eg /Users/akhil/work_dir/azure/AZURE2/azure2

mkdir -p "$HOST_OUTPUT" # ensure output directory exists

# --- Check for XQuartz ---
if ! command -v xquartz >/dev/null 2>&1; then
    echo " XQuartz not found. Install it from https://www.xquartz.org/"
    exit 1
fi

# --- Check for Docker ---
if ! command -v docker >/dev/null 2>&1; then
    echo "Docker not found. Install Docker Desktop from https://www.docker.com/products/docker-desktop"
    exit 1
fi

# --- Ensure Docker Desktop is running ---
if ! docker info >/dev/null 2>&1; then
    echo "🐳 Starting Docker Desktop..."
    open --background -a Docker
    echo -n "Waiting for Docker to start"
    until docker info >/dev/null 2>&1; do
        echo -n "."
        sleep 1
    done
    echo " Docker is ready."
fi

# --- Set LD_LIBRARY_PATH for ROOT inside the container ---
export LD_LIBRARY_PATH=/workspace/app/install/root/lib:/usr/lib:$LD_LIBRARY_PATH


# --- Ensure XQuartz allows local connections ---
xhost + 127.0.0.1 >/dev/null 2>&1 || true

# --- Check if the container already exists ---
if [ "$(docker ps -a -q -f name=$CONTAINER_NAME)" ]; then
    echo "🔁 Reattaching to existing container '$CONTAINER_NAME'..."
    docker start -ai $CONTAINER_NAME
else
    echo "🚀 Starting new container '$CONTAINER_NAME'..."
    docker run -it \
        --name $CONTAINER_NAME \
        --privileged \
        --hostname luna \
        -e SHELL=/bin/bash \
        -e DISPLAY=host.docker.internal:0 \
        -v "$HOME/.Xauthority:/root/.Xauthority" \
        -v "$APP_DIR:/workspace/azure2" \
        -v "$HOST_OUTPUT:/workspace/azure2/output" \
        azure2:latest bash -c '
            # Environment setup
            source /workspace/app/install/root/bin/thisroot.sh &&
            export QT_GRAPHICSSYSTEM="native" &&
            export QT_QUICK_BACKEND=software &&
            export LIBGL_ALWAYS_INDIRECT=1 &&

            # Ensure output folder exists
            mkdir -p /workspace/azure2/output &&

            # Launch AZURE2
            cd /workspace/azure2 &&
            AZURE2
        '
fi

echo "AZURE2 finished. Check your outputs at: $HOST_OUTPUT"
