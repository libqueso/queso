#!/bin/bash
set -eu
set -o pipefail

./test_jeffreys 
diff test_Regression/test_jeffreys_samples.m test_jeffreys_output/fp_p_seq.m
rm -r test_jeffreys_output
