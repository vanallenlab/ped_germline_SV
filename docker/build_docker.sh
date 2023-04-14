#!/usr/bin/env bash

# Study of Germline SVs in Pediatric Cancers
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Build project docker images

set -eu -o pipefail

TAG=$1

if [ -z $TAG ]; then
  echo "ERROR: Must provide desired image tag as a positional argument. Exiting."
  exit 1
fi

# Get various directories
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
EXEC_DIR=`pwd`
BUILD_DIR=`mktemp -d`

# Clone GATK-SV into build context
cd $BUILD_DIR && \
git clone git@github.com:broadinstitute/gatk-sv.git && \
cd $EXEC_DIR

# Copy ped SV repo into build context
cp -r $SCRIPT_DIR/../../ped_germline_SV $BUILD_DIR/

# Build base PedSV image
docker build \
  -f $SCRIPT_DIR/PedSV/Dockerfile \
  --progress plain \
  --tag vanallenlab/pedsv:$TAG \
  $BUILD_DIR

# Build PedSV-R image
docker build \
  -f $SCRIPT_DIR/PedSV-R/Dockerfile \
  --progress plain \
  --tag vanallenlab/pedsv-r:$TAG \
  $BUILD_DIR

# Push image & update latest
for image in pedsv pedsv-r; do
  docker push vanallenlab/${image}:$TAG
  docker tag vanallenlab/${image}:$TAG vanallenlab/${image}:latest
  docker push vanallenlab/${image}:latest
done
