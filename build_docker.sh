#!/bin/bash

time docker build -t my_alpine_build .
docker run --rm my_alpine_build