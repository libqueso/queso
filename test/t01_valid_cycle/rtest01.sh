#!/bin/bash
#----------------------------------------------------------
# PECOS Regression tests for QUESO
#
# 1. Validation cycle using TGA example.
#
# Originally: 5-19-09
#----------------------------------------------------------

set -eu
set -o pipefail

TOLERANCE="1e-10"	                   # solution diff tolerance (absolute)
SAVELOG=1		                   # Log model output?
COMMONDIR="$srcdir/common"

#----------------
# Initialization
#----------------

RUNDIR=`pwd`

. $COMMONDIR/verify.sh

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Regresion Test #1: Validation cycle with TGA example

VERIFY_DATE="04-06-2015"
SOLDIR="./output_test_TgaValidationCycle_gsl"
EXE="./TgaValidationCycle_gsl"   # executable name
SOLREFS="$srcdir/t01_valid_cycle/regression/$VERIFY_DATE"
TESTNAME='Test 1 (TGA Validation Cycle)'

# Setup desired input for the model

rm -f $SOLDIR/*.txt
rm -f $SOLDIR/*.m

verify_file_exists $EXE

# Remove output directory to nullify any previous results (and verify QUESO
# creates necessary parent directory).

rm -rf $SOLDIR

# Run the model

if [ $SAVELOG -eq 1 ];then
    ./$EXE
else
    ./$EXE >& /dev/null
fi

# Verify results

verify_file_exists $SOLDIR/file_val_ip_raw2.m
verify_file_exists $SOLREFS/file_val_ip_raw2.m

# Compare solutions from 4 output files.

for file in file_cal_ip_raw2.m file_val_ip_raw2.m file_cal_fp_qoi2.m file_val_fp_qoi2.m ; do
    if ! $RUNDIR/$COMMONDIR/compare.pl $SOLDIR/$file $SOLREFS/$file; then
      message_fail "$TESTNAME failed verification"
    fi
done

message_passed "$TESTNAME"
