#!/bin/bash

docker run -it \
       --privileged \
       --hostname luna \
       -e "DISPLAY=$DISPLAY" \
       -v /tmp/.X11-unix:/tmp/.X11-unix \
       -v "$XAUTH/.Xauthority:/root/.Xauthority" \
       -v "$PWD/azure2:/workspace/azure2" \
       azure2:latest bash -c 'cd /workspace/azure2 && \
       	      	                source /workspace/app/install/root/bin/thisroot.sh && \
				export QT_GRAPHICSSYSTEM="native" && \
       	      	       	        AZURE2'
