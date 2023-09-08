#!/usr/bin/env bash

# Study of Germline SVs in Pediatric Cancers
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Build project docker images

# Expects two positional arguments: image tag and (optional) comma-delimited list of images to build

set -eu -o pipefail

TAG=$1
if [ -z ${2:-} ]; then
  echo "No images specified; building all images by default"
  IMAGES="pedsv,pedsv-r"
else
  IMAGES=$2
fi

if [ -z $TAG ]; then
  echo "ERROR: Must provide desired image tag as a positional argument. Exiting."
  exit 1
fi

# Prune unused images before build (this is a common cause of failure)
docker image prune -f

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
if [ $( echo $IMAGES | sed 's/,/\n/g' | awk '{ if ($1=="pedsv") print }' | wc -l ) -gt 0 ]; then
  echo -e "\nPROGRESS: Now building PedSV image\n"
  docker build \
    -f $SCRIPT_DIR/PedSV/Dockerfile \
    --platform linux/x86_64 \
    --progress plain \
    --tag vanallenlab/pedsv:$TAG \
    $BUILD_DIR
fi

# Build PedSV-R image
if [ $( echo $IMAGES | sed 's/,/\n/g' | awk '{ if ($1=="pedsv-r") print }' | wc -l ) -gt 0 ]; then
  echo -e "\nPROGRESS: Now building PedSV-R image\n"
  docker build \
    -f $SCRIPT_DIR/PedSV-R/Dockerfile \
    --platform linux/x86_64 \
    --progress plain \
    --tag vanallenlab/pedsv-r:$TAG \
    $BUILD_DIR
fi

# Push image & update latest
for image in $( echo $IMAGES | sed 's/,/\n/g' ); do
  docker push vanallenlab/${image}:$TAG
  docker tag vanallenlab/${image}:$TAG vanallenlab/${image}:latest
  docker push vanallenlab/${image}:latest
done

# Clean up build images
docker image prune -f
