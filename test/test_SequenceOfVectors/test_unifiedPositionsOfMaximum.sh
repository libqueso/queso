#!/bin/bash
set -eu
set -o pipefail

PROG="./test_unifiedPositionsOfMaximum"

mpirun -np 2 $PROG
