#!/bin/bash
set -eu
set -o pipefail

PROG="./test_LlhdTargetOutput"

INPUT_1="${srcdir}/test_StatisticalInverseProblem/llhdout_input.txt"
INPUT_2="${srcdir}/test_StatisticalInverseProblem/targetout_input.txt"
INPUT_3="${srcdir}/test_StatisticalInverseProblem/both_input.txt"
INPUT_4="${srcdir}/test_StatisticalInverseProblem/neither_input.txt"

$PROG $INPUT_1
$PROG $INPUT_2
$PROG $INPUT_3
$PROG $INPUT_4

if [[ -f test_outputLlhd/ip_raw_chain_loglikelihood.m ]]
then
  echo "Found expected output file"
else
  exit 1
fi

if [[ -f test_outputTarget/ip_raw_chain_logtarget.m ]]
then
  echo "Found expected output file"
else
  exit 1
fi

if [[ -f test_outputBoth/ip_raw_chain_loglikelihood.m ]]
then
  echo "Found expected output file"
else
  exit 1
fi

if [[ -f test_outputBoth/ip_raw_chain_logtarget.m ]]
then
  echo "Found expected output file"
else
  exit 1
fi

if [[ -f test_outputNeither/ip_raw_chain_loglikelihood.m ]]
then
  echo "Found unexpected output file"
  exit 1
fi

if [[ -f tetst_outputNeither/ip_raw_chain_logtarget.m ]]
then
  echo "Found unexpected output file"
  exit 1
fi

rm -r test_outputLlhd
rm -r test_outputTarget
rm -r test_outputBoth
rm -r test_outputNeither
