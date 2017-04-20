#!/bin/bash
set -eu
set -o pipefail

TOLERANCE="1e-10"	 # solution diff tolerance (relative)
SAVELOG=1  # Log model output?
COMMONDIR="$srcdir/common"

RUNDIR=`pwd`

. $COMMONDIR/verify.sh

SOLDIR="./test_gpmsa_vector_output"
EXE="./test_gpmsa_vector"  # executable name
SOLREFS="${srcdir}/test_gpmsa"
INFILE="${srcdir}/test_gpmsa/gpmsa_vector_input.txt"
TESTNAME='Test Multivariate GPMSA'

rm -f $SOLDIR/*.m
rm -f $SOLDIR/*.txt

verify_file_exists $EXE
verify_file_exists $INFILE

rm -rf $SOLDIR

if [ $SAVELOG -eq 1 ];then
    ./$EXE $INFILE
else
    ./$EXE $INFILE >& /dev/null
fi

verify_file_exists $SOLDIR/ip_raw_chain.m
verify_file_exists $SOLREFS/test_gpmsa_vector_samples.m

$COMMONDIR/compare.pl $SOLDIR/ip_raw_chain.m $SOLREFS/test_gpmsa_vector_samples.m
exit $?
