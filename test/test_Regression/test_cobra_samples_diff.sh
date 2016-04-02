#!/bin/bash
set -eu
set -o pipefail

TOLERANCE="1e-10"	 # solution diff tolerance (relative)
SAVELOG=1  # Log model output?
COMMONDIR="$srcdir/common"

RUNDIR=`pwd`

. $COMMONDIR/verify.sh

SOLDIR="./test_gpmsa_cobra_output"
EXE="./test_gpmsa_cobra"  # executable name
SOLREFS="${srcdir}/test_Regression"
INFILE="${srcdir}/test_Regression/gpmsa_cobra_input.txt"
TESTNAME='Test Cobra GPMSA'

rm -f $SOLDIR/*.m
rm -f $SOLDIR/*.txt

verify_file_exists $EXE
verify_file_exists $INFILE

rm -rf ./test_gpmsa_cobra_output

if [ $SAVELOG -eq 1 ];then
    ./$EXE
else
    ./$EXE >& /dev/null
fi

verify_file_exists $SOLDIR/ip_raw_chain.m
verify_file_exists $SOLREFS/test_gpmsa_cobra_samples.m

$COMMONDIR/compare.pl $SOLDIR/ip_raw_chain.m $SOLREFS/test_gpmsa_cobra_samples.m
exit $?
