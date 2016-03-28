#!/bin/bash
set -eu
set -o pipefail

if (( @HAVE_MPI@ )); then
  mpirun -np 1 ../libtool --mode=execute ./test_intercomm0_gravity \
    test_intercomm0/gravity_1proc.txt

  mpirun -np 2 ../libtool --mode=execute ./test_intercomm0_gravity \
    test_intercomm0/gravity_2proc.txt

  for i in output_test_intercomm0_gravity_1/*.m; do
    j=`basename $i`;
    diff output_test_intercomm0_gravity_1/$j output_test_intercomm0_gravity_2/$j;
  done
else
  exit 77
fi
