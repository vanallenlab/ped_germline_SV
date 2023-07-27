#!/usr/bin/env bash

INVCF=$1
OUTVCF=$2

bcftools +fill-tags $INVCF -Oz -o $OUTVCF -- -t AC,AN,AF
