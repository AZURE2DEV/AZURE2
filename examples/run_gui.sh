#!/bin/bash

export XAUTH=$HOME/.Xauthority

docker run -it \
       --privileged \
       --hostname luna \
       -e DISPLAY=host.docker.internal:0 \
       -v /tmp/.X11-unix:/tmp/.X11-unix \
       -v "$XAUTH:/root/.Xauthority" \
       -v "$PWD/azure2:/workspace/azure2" \
       azure2:latest bash -c 'cd /workspace/azure2 && \
                                source /workspace/app/install/root/bin/thisroot.sh && \
                                export QT_GRAPHICSSYSTEM="native" && \
                                AZURE2'
