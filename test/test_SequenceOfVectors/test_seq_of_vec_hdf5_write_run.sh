#!/bin/bash
set -eu
set -o pipefail

../libtool --mode=execute ./test_seq_of_vec_hdf5_write
rm output_test_hdf5_sub0.h5
