#!/bin/bash
#----------------------------------------------------------
# PECOS Regression tests for QUESO
# 
# 1. Validation cycle using TGA example.
#
# Originally: 5-19-09
# 
# $Id: rtest.sh 2445 2009-04-21 08:19:13Z karl $
#----------------------------------------------------------

TOLERANCE="1e-10"	                   # solution diff tolerance (absolute)
SAVELOG=0		                   # Log model output?
COMMONDIR="$srcdir/../common"

#----------------
# Initialization
#----------------

RUNDIR=`pwd`

. $COMMONDIR/verify.sh

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Regresion Test #1: Validation cycle with TGA example 

VERIFY_DATE="07-30-2010"
SOLDIR="outputData"
EXE="./TgaValidationCycle_gsl"   # executable name
SOLREFS="$srcdir/regression/$VERIFY_DATE"
INFILE="$srcdir/tgaCycle.inp"
TESTNAME='Test 1 (TGA Validation Cycle)'

# Setup desired input for the model

rm -f $SOLDIR/*.txt
rm -f $SOLDIR/*.m

verify_file_exists $EXE
verify_file_exists $INFILE

# Remove output directory to nullify any previous results (and verify QUESO
# creates necessary parent directory).

rm -rf ./outputData

# Run the model

if [ $SAVELOG -eq 1 ];then
    ./$EXE $INFILE
else
    ./$EXE $INFILE >& /dev/null
fi

# Verify results

verify_file_exists $SOLDIR/file_val_ip_raw2.m
verify_file_exists $SOLREFS/file_val_ip_raw2.m

igot=1

# Compare solutions from 4 output files.

for file in file_cal_ip_raw2.m file_val_ip_raw2.m file_cal_fp_qoi2.m file_val_fp_qoi2.m ; do

    $RUNDIR/$COMMONDIR/compare.pl $SOLDIR/$file $SOLREFS/$file
    let igot="$igot * $?"
done

cd - >& /dev/null

if [ $igot -eq 0 ];then
    message_passed "$TESTNAME"
else
  message_fail "$TESTNAME failed verification"
fi

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------



