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
TOPDIR="./"		                   # relative to regression dir
SAVELOG=0		                   # Log model output?

#----------------
# Initialization
#----------------

RUNDIR=`pwd`

. ./t01_valid_cycle/verify.sh

cd $TOPDIR
verify_file_exists $EXE

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/h2/prudenci/Installations/Boost_1_35_0/lib/

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Regresion Test #1: Validation cycle with TGA example 

VERIFY_DATE="05-19-2009"
TEST_DIR="t01_valid_cycle/valid_cycle"
SOLDIR="outputData"
EXE="./TgaValidationCycle_gsl"   # executable name
SOLREFS="regression/$VERIFY_DATE"
INFILE="tgaCycle.inp"
TESTNAME='Test 1 (TGA Validation Cycle)'

# Setup desired input for the model

cd $TEST_DIR

rm -f $SOLDIR/*.m
verify_file_exists $INFILE

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

    $RUNDIR/t01_valid_cycle/compare.pl $SOLDIR/$file $SOLREFS/$file
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



