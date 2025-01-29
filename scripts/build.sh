#!/bin/bash

docker build --ulimit nofile=1024 -t azure2 .
