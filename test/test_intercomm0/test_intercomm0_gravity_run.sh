#!/bin/bash
set -eu
set -o pipefail

if grep "QUESO_HAVE_MPI 1" ../config_queso.h 2>&1 >/dev/null; then
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
