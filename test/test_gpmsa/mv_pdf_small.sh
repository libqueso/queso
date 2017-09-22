#!/bin/bash
set -eu
set -o pipefail

./test_gpmsa_scalar_pdf_small ${srcdir}/test_gpmsa/gpmsa_mv_pdf_small_input.txt
