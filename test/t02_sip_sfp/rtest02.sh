#!/bin/bash
#----------------------------------------------------------
# PECOS Regression tests for QUESO
#
# 1. Validation cycle using TGA example.
#
# Originally: 5-19-09
#----------------------------------------------------------

TOLERANCE="1e-10"	                   # solution diff tolerance (relative)
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

#VERIFY_DATE="07-12-2009"
#VERIFY_DATE="09-06-2009"
#VERIFY_DATE="07-09-2012"
#VERIFY_DATE="09-03-2012"
#VERIFY_DATE="09-04-2012"
# VERIFY_DATE="06-11-2014"
# VERIFY_DATE="09-15-2014"
VERIFY_DATE="05-01-2015"
SOLDIR="./output_test_SipSfpExample_gsl"
EXE="./SipSfpExample_gsl"   # executable name
SOLREFS="$srcdir/t02_sip_sfp/regression/$VERIFY_DATE"
INFILE="$srcdir/t02_sip_sfp/example.inp"
TESTNAME='Test 2 (Sip + Sfp example)'

# Setup desired input for the model

rm -f $SOLDIR/*.txt
rm -f $SOLDIR/*.m

verify_file_exists $EXE
verify_file_exists $INFILE

# Remove output directory to nullify any previous results (and verify QUESO
# creates necessary parent directory).

rm -rf $SOLDIR

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
    grep "e-[0-9]\|e+[0-9]" $SOLREFS/$file | grep -v sec | grep -v Arch | grep -v "Build Host" > $SOLDIR/nada0
    grep "e-[0-9]\|e+[0-9]" $SOLDIR/$file | grep -v sec | grep -v Arch | grep -v "Build Host" > $SOLDIR/nada1
    diff $SOLDIR/nada0 $SOLDIR/nada1
    let igot="$igot + $?"
done

cd - >& /dev/null

if [ $igot -eq 0 ];then
    message_passed "$TESTNAME"
else
  message_fail "$TESTNAME failed verification"
fi
