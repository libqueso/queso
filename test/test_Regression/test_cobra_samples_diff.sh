#!/bin/bash
set -eu
set -o pipefail

./test_gpmsa_cobra
diff test_Regression/test_gpmsa_cobra_samples.m test_gpmsa_cobra_output/ip_raw_chain_sub0.m
