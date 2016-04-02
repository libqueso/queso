#!/bin/bash
set -eu
set -o pipefail

if [ -d test_jeffreys_output ]; then
  rm -r test_jeffreys_output
fi
./test_jeffreys 
diff ${srcdir}/test_Regression/test_jeffreys_samples.m test_jeffreys_output/fp_p_seq.m
rm -r test_jeffreys_output
