#!/bin/bash
set -eu
set -o pipefail

./test_gpmsa_scalar_pdf_large ${srcdir}/test_gpmsa/gpmsa_scalar_pdf_large_input.txt
