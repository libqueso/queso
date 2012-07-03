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
COMMONDIR="../common"

#----------------
# Initialization
#----------------

RUNDIR=`pwd`

. $COMMONDIR/verify.sh

cd $TOPDIR
verify_file_exists $EXE

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
# Regresion Test #1: Validation cycle with TGA example 

#VERIFY_DATE="07-12-2009"
#VERIFY_DATE="09-22-2010"
#VERIFY_DATE="12-09-2010"
VERIFY_DATE="07-03-2012"
TEST_DIR="./"
SOLDIR="outputData"
EXE="./BimodalExample_gsl"   # executable name
SOLREFS="regression/$VERIFY_DATE"
INFILE="example_1chain.inp"
TESTNAME='Test 4 (bimodal example)'

# Setup desired input for the model

cd $TEST_DIR

rm -f $SOLDIR/*.txt
rm -f $SOLDIR/*.m
verify_file_exists $INFILE

# Run the model

if [ $SAVELOG -eq 1 ];then
    ./$EXE $INFILE
else
    ./$EXE $INFILE >& /dev/null
fi

# Verify results

verify_file_exists $SOLDIR/display_sub0.txt
verify_file_exists $SOLREFS/display_sub0.txt

igot=0

# Compare outputs from 1 output file.

for file in display_sub0.txt ; do

#    pwd
#    echo $SOLREFS/$file
    grep "e-\|e+" $SOLREFS/$file | grep -v sec | grep -v Arch > $SOLDIR/nada0
    grep "e-\|e+" $SOLDIR/$file | grep -v sec | grep -v Arch > $SOLDIR/nada1
    diff $SOLDIR/nada0 $SOLDIR/nada1
    let igot="$igot + $?"

done

cd - >& /dev/null

if [ $igot -eq 0 ];then
    message_passed "$TESTNAME"
else
  message_fail "$TESTNAME failed verification"
fi

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------



