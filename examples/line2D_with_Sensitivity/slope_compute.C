/* This file performs the statistical inference for the slope
 * of the straight line: y = mx. In addition to 'compute_slope.h',
 * it relies on 'slope_likelihood.h', corresponding to 'slope_likelihood.C'
 * file.
 */

#include <cmath>
#include <sys/time.h>
#include <stdlib.h>

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>

#include <slope_compute.h>
#include <slope_likelihood.h>
#include <slope_qoi.h>
#include <sensitivity_m.h>
#include <sensitivity_c.h>
#include <sensitivity_mc.h>

void infer_slope(const QUESO::FullEnvironment & env) {


// Statistical Inverse Problem: Compute posterior pdf for slope 'm' and y-intercept c

// Step 1: Instantiate the parameter space

	QUESO::VectorSpace<> paramSpace(env, "param_", 2, NULL); // 2 since we have a 2D problem
	
// Step 2: Parameter domain

	QUESO::GslVector paramMinValues(paramSpace.zeroVector());
	QUESO::GslVector paramMaxValues(paramSpace.zeroVector());

	paramMinValues[0] = 2.;
	paramMaxValues[0] = 5.;

	paramMinValues[1] = 3.;
	paramMaxValues[1] = 7.;
	
	QUESO::BoxSubset<> paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);

// Step 3: Instantiate likelihood

	Likelihood<> lhood("like_", paramDomain);

// Step 4: Define the prior RV

	QUESO::UniformVectorRV<> priorRv("prior_", paramDomain);

// Step 5: Instantiate the inverse problem

	QUESO::GenericVectorRV<> postRv("post_", paramSpace);
	QUESO::StatisticalInverseProblem<> ip("", NULL, priorRv, lhood, postRv);

// Step 6: Solve the inverse problem

// Randomly sample for the initial state?
	QUESO::GslVector paramInitials(paramSpace.zeroVector());
	priorRv.realizer().realization(paramInitials);

// Initialize the Cov matrix:
	QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
	proposalCovMatrix(0,0) = std::pow(std::abs(paramInitials[0]) / 20.0, 2.0);
	proposalCovMatrix(1,1) = std::pow(std::abs(paramInitials[1]) / 20.0, 2.0);

	ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

// Using the posterior pdfs for m and c, compute 'y' at a given 'x'
 
// Step 1: Instantiate the qoi space

	QUESO::VectorSpace<> qoiSpace(env, "qoi_", 1, NULL);

// Step 2: Instantiate the parameter domain

// Not necessary here because the posterior from SIP is used as the RV for SFP

// Step 3: Instantiate the qoi object to be used by QUESO

	Qoi<> qoi("qoi_", paramDomain, qoiSpace);

// Step 4: Define the input RV

// Not required because we use the posterior as RV

// Step 5: Instantiate the forward problem

	QUESO::GenericVectorRV<> qoiRv("qoi_", qoiSpace);
	
	QUESO::StatisticalForwardProblem<> fp("", NULL, postRv, qoi, qoiRv);

// Step 6: Solve the forward problem

	std::cout << "Solving the SFP with Monte Carlo" << std::endl << std::endl;
	fp.solveWithMonteCarlo(NULL);

	system("mv outputData/sfp_lineSlope_qoi_seq.txt outputData/sfp_lineSlope_qoi_seq_post.txt");

// SENSITIVITY ANALYSIS

// For m

	Qoi_m<> qoi_m("qoi_", paramDomain, qoiSpace);
		
// Step 4: Define the input RV

// Not required because we use the prior as RV for sensitivity analysis

// Step 5: Instantiate the forward problem

	QUESO::StatisticalForwardProblem<> fp_m("", NULL, priorRv, qoi_m, qoiRv);

// Step 6: Solve the forward problem

	fp_m.solveWithMonteCarlo(NULL);

	system("mv outputData/sfp_lineSlope_qoi_seq.txt outputData/sense_m.txt");

// For c

	Qoi_c<> qoi_c("qoi_", paramDomain, qoiSpace);
		
// Step 4: Define the input RV

// Not required because we use the prior as RV for sensitivity analysis

// Step 5: Instantiate the forward problem

	QUESO::StatisticalForwardProblem<> fp_c("", NULL, priorRv, qoi_c, qoiRv);

// Step 6: Solve the forward problem

	fp_c.solveWithMonteCarlo(NULL);

	system("mv outputData/sfp_lineSlope_qoi_seq.txt outputData/sense_c.txt");

// For both, m and c

	Qoi_mc<> qoi_mc("qoi_", paramDomain, qoiSpace);
		
// Step 4: Define the input RV

// Not required because we use the prior as RV for sensitivity analysis

// Step 5: Instantiate the forward problem

	QUESO::StatisticalForwardProblem<> fp_mc("", NULL, priorRv, qoi_mc, qoiRv);

// Step 6: Solve the forward problem

	fp_mc.solveWithMonteCarlo(NULL);

	system("mv outputData/sfp_lineSlope_qoi_seq.txt outputData/sense_mc.txt");
}
