#!/bin/sh
set -e

(cd gp/scalar && ./gpmsa_scalar gpmsa_input.txt)
(cd bimodal && ./bimodal_gsl bimodal_1chain.inp)
(cd gravity && ./gravity_gsl gravity_inv_fwd.inp)
(cd hysteretic && ./hysteretic_gsl example.inp)
(cd infinite_dim && ./inverse_problem inverse_options.in)
(cd infinite_dim &&
  mpirun -np 2 ./parallel_inverse_problem parallel_inverse_options.in)

# Only works with --enable-ann
(cd infoTheoryProblem && ./exInfoTheory_gsl infTh.inp || echo "Skipping example")

(cd interpolation_surrogate && ./4d_interp input.in)
(cd simpleStatisticalForwardProblem &&
  ./exSimpleStatisticalForwardProblem_gsl simple_sfp_example.inp)
(cd simpleStatisticalInverseProblem &&
  ./exSimpleStatisticalInverseProblem_gsl example.inp)
(cd validationCycle && ./exTgaValidationCycle_gsl tgaCycle.inp)
(cd validationCycle && ./exTgaValidationCycle_gsl tgaCycle_2009_03_30.inp)
(cd validationCycle2 && ./tga2_gsl tgaCycle.inp)
