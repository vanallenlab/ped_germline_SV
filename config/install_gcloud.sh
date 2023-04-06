#!/usr/bin/env bash

# Install Google Cloud SDK on Ubuntu

apt-get -qqy clean && \
apt-get -qqy update --fix-missing && \
apt-get -qqy dist-upgrade && \
apt-get -qqy install lsb-release apt-transport-https ca-certificates gnupg && \
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
| tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
apt-get -qqy update && \
apt-get -qqy install --no-install-recommends google-cloud-sdk && \
gcloud config set core/disable_usage_reporting true && \
gcloud config set component_manager/disable_update_check true && \
gcloud config set metrics/environment github_docker_image && \
apt-get -qqy clean && \
rm -rf /tmp/* \
       /var/tmp/* \
       /var/cache/apt/* \
       /var/lib/apt/lists/* \
       /usr/share/man/?? \
       /usr/share/man/??_* && \
gcloud --help
