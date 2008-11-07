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

#ifndef __UQ_BURGERS_NO_MODEL_UNCERTAINTY_H__
#define __UQ_BURGERS_NO_MODEL_UNCERTAINTY_H__

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

// The (user defined) data class for the data needed by the 
// (user defined) likelihood routine
template<class P_V, class P_M>
struct 
likelihoodRoutine_DataClass
{
  likelihoodRoutine_DataClass(const uqEnvironmentClass& env,
                              const char* fileName);
 ~likelihoodRoutine_DataClass();

  // Experimental data
  int nDataPoints; // number of data points
  double *dataLocations; // x locations where state was measured
  double *dataValues; // measured values of averaged state at dataLocations

  // For Burgers' solver
  quadBasis *pQB; // quad points, weights, and basis evaluation
  gsl_vector *U; // solution (re-used after each solve at IC for next solve to hopefully decrease Newton iterations)
};

// likelihoodRoutine_DataClass constructor: Read data, allocate/initialize memory for Burgers' solver
template<class P_V, class P_M>
likelihoodRoutine_DataClass<P_V,P_M>::likelihoodRoutine_DataClass(const uqEnvironmentClass &env, const char *fileName)
  :
  nDataPoints(0),
  dataLocations(0),
  dataValues(0),
  pQB(0),
  U(0)
{
  std::cout << "Calling likelihoodRoutine_DataClass constructor...";
  
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
  int nDataPoints = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->nDataPoints;
  const double *xx = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->dataLocations;
  const double *ue = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->dataValues;

  // Solver data
  quadBasis *pQB = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->pQB;
  gsl_vector *U = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->U;

  // Model parameter
  const double kappa = paramValues[0];

  // Model data
  double um[nDataPoints];

  
  // Evaluate model
  int ierr;
  ierr = solveForStateAtXLocations(xx, nDataPoints, 1e-1, kappa, U, pQB, um);
  if( ierr != 0 ){
    printf("WARNING: Burgers solver was not successful!\n");
    printf("         Results are probably meaningless.\n");
    fflush(stdout);
  }

  // Compute misfit vector
  double misfit2[nDataPoints];
  for( int ii=0; ii<nDataPoints; ii++ ) misfit2[ii] = (ue[ii] - um[ii])*(ue[ii] - um[ii]);

  double sum = 0.0;
  for( int ii=0; ii<nDataPoints; ii++ ) sum += misfit2[ii];

  double var = 1e-8;
  double resultValue = sum/var;

  return resultValue;
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

  int iRC;
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

  likelihoodRoutine_DataClass<P_V,P_M> calLikelihoodRoutine_Data(env, "calibration_data.dat");

  cycle.setCalIP(calPriorRv,
		 likelihoodRoutine<P_V,P_M>,
		 (void *) &calLikelihoodRoutine_Data,
		 true); // the likelihood routine computes [-2.*ln(Likelihood)]


  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv' (using Markov Chain)
  P_M* calProposalCovMatrix = 
    cycle.calIP().postRv().imageSpace().newGaussianMatrix(cycle.calIP().priorRv().pdf().domainVarianceValues(),
							  calInitialValues);

  cycle.calIP().solveWithBayesMarkovChain(calInitialValues,
                                          *calProposalCovMatrix,
                                          NULL); // use default kernel from library
  delete calProposalCovMatrix;


  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqBurgers_No_Model_Uncertainty' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_BURGERS_NO_MODEL_UNCERTAINTY_H__
