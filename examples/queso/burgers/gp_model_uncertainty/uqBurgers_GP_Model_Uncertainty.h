/* uq/examples/queso/tga/uqTgaEx4.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_BURGERS_GP_MODEL_UNCERTAINTY_H__
#define __UQ_BURGERS_GP_MODEL_UNCERTAINTY_H__

#define PI 3.14159265358979323846

// queso
#include <uqDefines.h>
#include <uqEnvironment.h>
#include <uqVectorSpace.h>
#include <uqVectorRV.h>
#include <uqAsciiTable.h>
#include <uqValidationCycle.h>

// system
#include <sys/time.h>

// gsl
#include <gsl/gsl_vector.h>

// burgers
#include <tools.h>
#include <quesoFunctions.h>

//********************************************************
// Likelihood function object for both inverse problems of the validation cycle.
// A likelihood function object is provided by user and is called by the UQ library.
// This likelihood function object consists of data and routine.
//********************************************************

// Declaration of function to compute likelihood from misfit
// (defined in uqBurgers_GP_Model_Uncertainty.C)
void 
negTwoLogLikelihood(const int N, const double *cholK, const double *misfit, 
		    double *result);

// Declaration of function to compute GP covariance matrix for calibration
// (defined in uqBurgers_GP_Model_Uncertainty.C)
int
calibrationCovarianceChol(const int N, const double *xLoc, const double sig2, const double ellx, 
			  double *cholK);


// The (user defined) data class for the data needed by the 
// (user defined) likelihood routine
template<class P_V, class P_M>
struct 
likelihoodRoutine_DataClass
{
  likelihoodRoutine_DataClass(const uqEnvironmentClass& env,
			      double Re,
                              const char* fileName);
 ~likelihoodRoutine_DataClass();

  // Experimental data
  int nDataPoints; // number of data points
  double Re;
  double *dataLocations; // x locations where state was measured
  double *dataValues; // measured values of averaged state at dataLocations

  // Gaussian process model uncertainty
  double *mean;  // GP mean evaluated at data points
  double *cholK; // GP covariance evaluated at data points

  // For Burgers' solver
  quadBasis *pQB; // quad points, weights, and basis evaluation
  gsl_vector *U; // solution (re-used after each solve at IC for next solve to hopefully decrease Newton iterations)
};

// likelihoodRoutine_DataClass constructor: Read data, allocate/initialize memory for Burgers' solver
template<class P_V, class P_M>
likelihoodRoutine_DataClass<P_V,P_M>::likelihoodRoutine_DataClass(const uqEnvironmentClass &env, 
								  double Re, const char *fileName)
  :
  nDataPoints(0),
  dataLocations(0),
  dataValues(0),
  pQB(0),
  U(0)
{
  std::cout << "Calling likelihoodRoutine_DataClass constructor...";

  // set Re
  this->Re = Re;
  
  // open file
  FILE *fp = fopen(fileName, "r");
  UQ_FATAL_TEST_MACRO((!fp), env.rank(),
		      "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		      "could not open input file");
  
  // get number of data points
  UQ_FATAL_TEST_MACRO( (fscanf(fp, "%d", &nDataPoints)!=1), env.rank(),
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		       "could not read number of data points");

  // allocate vectors
  dataLocations = new double[nDataPoints];
  dataValues = new double[nDataPoints];

  // read data
  int numObservations=0;
  double tmpX, tmpU;
  while (fscanf(fp,"%lf %lf",&tmpX,&tmpU) != EOF) {
    UQ_FATAL_TEST_MACRO((numObservations >= nDataPoints), env.rank(),
			"uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
			"failed while reading data");
    dataLocations[numObservations] = tmpX;
    dataValues[numObservations]    = tmpU;
    numObservations++;
  }

  // close file
  fclose(fp);

  // allocate and initialize mean vector and covariance matrix for Gaussian process
  mean = new double[nDataPoints];
  for( int ii=0; ii<nDataPoints; ii++ ) mean[ii] = 0.0;

  cholK = new double[nDataPoints*nDataPoints];
  for( int ii=0; ii<nDataPoints; ii++ ){
    for( int jj=0; jj<nDataPoints; jj++ ){
      cholK[nDataPoints*ii+jj] = 0.0;
    }
    cholK[nDataPoints*ii+ii] = 1.0;
  }

  // allocate and initialize quadrature points
  UQ_FATAL_TEST_MACRO( (evaluateQuadratureAndBasisForResidual(100, &pQB)!=0), env.rank(),
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		       "failed while allocating/computing quadrature points" );
		       
  // allocate and initialize solution (to zero) 
  U = gsl_vector_calloc(100);

  std::cout << "success.\n";
} // end constructor


// likelihoodRoutine_DataClass destructor: Frees memory allocated by constructor
template<class P_V, class P_M>
likelihoodRoutine_DataClass<P_V,P_M>::~likelihoodRoutine_DataClass()
{
  std::cout << "Calling likelihoodRoutine_DataClass destructor...";
  delete[] dataLocations;
  delete[] dataValues;
  delete[] mean;
  delete[] cholK;
  freeQuadratureAndBasisForResidual(pQB);
  gsl_vector_free(U);
  std::cout << "success.\n";
} // end destructor


// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
likelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  //std::cout << "Calling likelihoodRoutine...\n";

  // Experimental data
  const int nDataPoints = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->nDataPoints;
  const double Re = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->Re;
  const double *xx = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->dataLocations;
  const double *ue = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->dataValues;

  // Gaussian process model uncertainty
  const double *mean = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->mean;
  const double *cholK = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->cholK;

  // Solver data
  quadBasis *pQB = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->pQB;
  gsl_vector *U = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->U;

  // Model parameter
  const double kappa = paramValues[0];

  // Model data
  double um[nDataPoints];

  
  // Evaluate model
  int ierr = solveForStateAtXLocations(xx, nDataPoints, 1.0/Re, kappa, U, pQB, um);
  UQ_FATAL_TEST_MACRO( (ierr!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine)",
		       "failed solving Burgers' eqn" );

  // Compute misfit vector
  double misfit[nDataPoints];
  for( int ii=0; ii<nDataPoints; ii++ ) misfit[ii] = (ue[ii] - um[ii] - mean[ii]);

  double resultValue;
  negTwoLogLikelihood(nDataPoints, cholK, misfit, &resultValue);

  return resultValue;
}


//********************************************************
// QoI function object for the first validation problem stage (with prefix "s1_").
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data type for the data needed by the (user defined) qoi routine
struct
qoiRoutine_DataClass
{
  qoiRoutine_DataClass();
  ~qoiRoutine_DataClass();

  // For the Burgers solver
  quadBasis *pQB; // quad points, weights, and basis evaluation
};

// constructor
qoiRoutine_DataClass::qoiRoutine_DataClass()
{
  cout << "Calling qoiRoutine_DataClass constructor...";
  // allocate and initialize quadrature points
  UQ_FATAL_TEST_MACRO( (evaluateQuadratureAndBasisForResidual(100, &pQB)!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		       "failed while allocating/computing quadrature points" );
  cout << "success.\n";
} // end constructor

// destructor
qoiRoutine_DataClass::~qoiRoutine_DataClass()
{
  cout << "Calling qoiRoutine_DataClass destructor...";
  freeQuadratureAndBasisForResidual(pQB);
  cout << "success.\n";
} // end destructor

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void qoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{

  int ierr;
  const double kappa = paramValues[0];
  double u_x1;
  quadBasis *pQB = ((qoiRoutine_DataClass *) functionDataPtr)->pQB;

  UQ_FATAL_TEST_MACRO( (computeGradientAtOne(kappa, pQB, &u_x1)!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (qoiRoutine)",
		       "failed while computing QoI.");
  qoiValues[0] = u_x1;
  return;
}


//********************************************************
// The driving routine "uqAppl()": called by main()
// Stage   I: the 'calibration stage'
// Stage  II: the 'validation stage'
// Stage III: the 'comparison stage'
//********************************************************
template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqAppl(const uqEnvironmentClass& env)
{
  if (env.rank() == 0) {
    std::cout << "Beginning run of 'uqBurgers_No_Model_Uncertainty' example\n"
              << std::endl;
  }

  UQ_FATAL_TEST_MACRO(env.isThereInputFile() == false,
                      env.rank(),
                      "uqAppl()",
                      "input file must be specified in command line, after the '-i' option");

  int iRC, ierr;
  struct timeval timevalRef;
  struct timeval timevalNow;

  //******************************************************
  // Burgers' example: instantiation
  //******************************************************

  // Read Ascii file with model parameter information
  uqAsciiTableClass<P_V,P_M> paramTable(env,
					1,    // # of rows
					3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
					NULL, // All extra columns are of 'double' type
					"param.tab");


  const EpetraExt::DistArray<std::string>& paramNames = paramTable.stringColumn(0);
  P_V paramMinValues    (paramTable.doubleColumn(1));
  P_V paramMaxValues    (paramTable.doubleColumn(2));
  P_V calInitialValues  (paramTable.doubleColumn(3));

  uqVectorSpaceClass<P_V,P_M> paramSpace(env,
                                         "param_", // Extra prefix before the default "space_" prefix
                                         paramTable.numRows(),
                                         &paramNames);

  // Read Ascii file with QoI information.
  uqAsciiTableClass<P_V,P_M> qoiTable(env,
                                      1,    // # of rows
                                      0,    // # of cols after 'parameter name': none
                                      NULL, // All extra columns are of 'double' type
                                      "qoi.tab");

  const EpetraExt::DistArray<std::string>& qoiNames = qoiTable.stringColumn(0);

  uqVectorSpaceClass<Q_V,Q_M> qoiSpace(env,
                                       "qoi_", // Extra prefix before the default "space_" prefix
                                       qoiTable.numRows(),
                                       &qoiNames);

  // Instantiate the validation cycle
  uqValidationCycleClass<P_V,P_M,Q_V,Q_M> cycle(env,
                                                "", // No extra prefix
                                                paramSpace,
                                                qoiSpace);

  //******************************************************
  // Burgers' example: calibration phase
  //******************************************************
  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning 'calibration stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Set up inverse problem (i.e. prior and likelihood)
  uqUniformVectorRVClass<P_V,P_M> calPriorRv("cal_prior_", // Extra prefix before the default "rv_" prefix
                                             paramSpace,
                                             paramMinValues,
                                             paramMaxValues);

  likelihoodRoutine_DataClass<P_V,P_M> calLikelihoodRoutine_Data(env, 10.0, "calibration_data.dat");

  // Set GP covariance for calibration phase (from prior)
  ierr = calibrationCovarianceChol(calLikelihoodRoutine_Data.nDataPoints,
				   calLikelihoodRoutine_Data.dataLocations,
				   1e-4, 0.1,
				   calLikelihoodRoutine_Data.cholK);
  UQ_FATAL_TEST_MACRO( (ierr!=0), env.rank(),
		       "uqAppl(), in uqBurgers_GP_Model_Uncertainty",
		       "failed while computing calibration covariance matrix" );


  cycle.setCalIP(calPriorRv, likelihoodRoutine<P_V,P_M>, (void *) &calLikelihoodRoutine_Data,
		 true); // the likelihood routine computes [-2.*ln(Likelihood)]


  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv' (using Markov Chain)
  P_M* calProposalCovMatrix = 
    cycle.calIP().postRv().imageSpace().newGaussianMatrix(cycle.calIP().priorRv().pdf().domainVarianceValues(),
							  calInitialValues);

  cycle.calIP().solveWithBayesMarkovChain(calInitialValues,
                                          *calProposalCovMatrix,
                                          NULL); // use default kernel from library
  delete calProposalCovMatrix;



  // done with calibration phase
  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending 'calibration stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'calibration stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  //******************************************************
  // Burgers' example: validation phase
  //******************************************************
  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning 'validation stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

/*   // Set up inverse problem */
/*   likelihoodRoutine_DataClass<P_V,P_M> valLikelihoodRoutine_Data(env, 100.0, "validation_data.dat"); */

/*   cycle.setValIP(likelihoodRoutine<P_V,P_M>, (void *) &valLikelihoodRoutine_Data, */
/*                  true); // the likelihood routine computes [-2.*ln(Likelihood)] */


/*   // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv' */
/*   // Use 'realizer()' because the posterior rv was computed with Markov Chain */
/*   // Use these calibration mean values as the initial values */
/*   P_M* valProposalCovMatrix =  */
/*     cycle.calIP().postRv().imageSpace().newGaussianMatrix(cycle.calIP().postRv().realizer().imageVarianceValues(), */
/* 							  cycle.calIP().postRv().realizer().imageExpectedValues());  */

/*   cycle.valIP().solveWithBayesMarkovChain(cycle.calIP().postRv().realizer().imageExpectedValues(), */
/*                                           *valProposalCovMatrix, */
/*                                           NULL); // use default kernel from library */
/*   delete valProposalCovMatrix; */


  // done with validation phase
  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending 'calibration stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'calibration stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }


  //******************************************************
  // Burgers' example: comparison phase
  //******************************************************

  // done with everything
  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqBurgers_No_Model_Uncertainty' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_BURGERS_GP_MODEL_UNCERTAINTY_H__
