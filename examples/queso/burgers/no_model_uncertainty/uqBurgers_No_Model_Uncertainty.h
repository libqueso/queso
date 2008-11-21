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
#include <quesoInterface.h>

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
  likelihoodRoutine_DataClass(const uqBaseEnvironmentClass& env,
			      double Re,
                              const char* fileName,
			      burgersQuesoInterface& interface);
 ~likelihoodRoutine_DataClass();

  // burgers/queso interface class
  burgersQuesoInterface& burgersInterface;

  // Experimental data
  int nDataPoints; // number of data points
  double Re;
  double *dataLocations; // x locations where state was measured
  double *dataValues; // measured values of averaged state at dataLocations

};

// likelihoodRoutine_DataClass constructor: Read data, allocate/initialize memory for Burgers' solver
template<class P_V, class P_M>
likelihoodRoutine_DataClass<P_V,P_M>::likelihoodRoutine_DataClass(const uqBaseEnvironmentClass &env, 
								  double Re, const char *fileName,
								  burgersQuesoInterface& interface)
  :
  burgersInterface(interface),
  nDataPoints(0),
  dataLocations(0),
  dataValues(0)
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

  std::cout << "success.\n";
} // end constructor


// likelihoodRoutine_DataClass destructor: Frees memory allocated by constructor
template<class P_V, class P_M>
likelihoodRoutine_DataClass<P_V,P_M>::~likelihoodRoutine_DataClass()
{
  std::cout << "Calling likelihoodRoutine_DataClass destructor...";
  delete[] dataLocations;
  delete[] dataValues;
  std::cout << "success.\n";
} // end destructor


// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
likelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  // Experimental data
  const int nDataPoints = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->nDataPoints;
  const double Re = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->Re;
  const double *ue = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->dataValues;
  burgersQuesoInterface& burgersInterface = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->burgersInterface;

  // Model parameter
  double kappa = paramValues[0];

  // Model data
  double um[nDataPoints];

  // Evaluate model
  int ierr = burgersInterface.solveForStateAtDataLocations(1.0/Re, kappa, um);
  UQ_FATAL_TEST_MACRO( (ierr!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine)",
		       "failed solving Burgers' eqn" );

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
// QoI function object for the first validation problem stage (with prefix "s1_").
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data type for the data needed by the (user defined) qoi routine
struct
qoiRoutine_DataClass
{
  qoiRoutine_DataClass(burgersQuesoInterface& interface);
  ~qoiRoutine_DataClass();

  double qoiRe;

  // burgers/queso interface class
  burgersQuesoInterface& burgersInterface;

};

// constructor
qoiRoutine_DataClass::qoiRoutine_DataClass(burgersQuesoInterface& interface)
  :
  burgersInterface(interface)
{
  cout << "Calling qoiRoutine_DataClass constructor...";
  qoiRe = 200.0;
  cout << "success.\n";
} // end constructor

// destructor
qoiRoutine_DataClass::~qoiRoutine_DataClass()
{
  cout << "Calling qoiRoutine_DataClass destructor...";
  cout << "success.\n";
} // end destructor

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void qoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{

  int ierr;
  double kappa = paramValues[0];
  double u_x1;
  double qoiRe = ((qoiRoutine_DataClass *) functionDataPtr)->qoiRe;
  burgersQuesoInterface& burgersInterface = ((qoiRoutine_DataClass *) functionDataPtr)->burgersInterface;

  ierr = burgersInterface.solveForViscFluxAtOne(1.0/qoiRe, kappa, &u_x1);
  UQ_FATAL_TEST_MACRO( (ierr!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (qoiRoutine)",
		       "failed solving Burgers' eqn" );

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
uqAppl(const uqBaseEnvironmentClass& env)
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

  uqBoxSubsetClass<P_V,P_M> paramDomain("param_",
                                        paramSpace,
                                        paramMinValues,
                                        paramMaxValues);

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
                                             paramDomain);

  burgersQuesoInterface cal_burgersInterface(100, 10.0, "calibration_data.dat");
  likelihoodRoutine_DataClass<P_V,P_M> calLikelihoodRoutine_Data(env, 10.0, "calibration_data.dat", cal_burgersInterface);


  uqGenericScalarFunctionClass<P_V,P_M> calLikelihoodFunctionObj("cal_like_",
                                                                 paramDomain,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 NULL,
                                                                 NULL,
                                                                 (void *) &calLikelihoodRoutine_Data,
                                                                 true);
  
  cycle.setCalIP(calPriorRv, calLikelihoodFunctionObj);
  
  P_M* calProposalCovMatrix = 
    cycle.calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(cycle.calIP().priorRv().pdf().domainVarVector(), 
								      calInitialValues);

  cycle.calIP().solveWithBayesMarkovChain(calInitialValues, calProposalCovMatrix);
  delete calProposalCovMatrix;

  // Deal with forward problem
  qoiRoutine_DataClass calQoiRoutine_Data(cal_burgersInterface);

  cycle.setCalFP(qoiRoutine<P_V,P_M,Q_V,Q_M>, (void *) &calQoiRoutine_Data);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  cycle.calFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

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

  // Set up inverse problem
  burgersQuesoInterface val_burgersInterface(100, 100.0, "validation_data.dat");
  likelihoodRoutine_DataClass<P_V,P_M> valLikelihoodRoutine_Data(env, 100.0, "validation_data.dat", val_burgersInterface);

  uqGenericScalarFunctionClass<P_V,P_M> valLikelihoodFunctionObj("val_like_",
                                                                 paramDomain,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 NULL,
                                                                 NULL,
                                                                 (void *) &valLikelihoodRoutine_Data,
                                                                 true); // the routine computes [-2.*ln(function)]

  cycle.setValIP(valLikelihoodFunctionObj);


   P_M* valProposalCovMatrix = 
     cycle.calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(cycle.calIP().postRv().realizer().imageVarVector(), cycle.calIP().postRv().realizer().imageExpVector()); 

  cycle.valIP().solveWithBayesMarkovChain(cycle.calIP().postRv().realizer().imageExpVector(),
                                          valProposalCovMatrix);
  delete valProposalCovMatrix;


  // Deal with forward problem
  qoiRoutine_DataClass valQoiRoutine_Data(val_burgersInterface);

  cycle.setValFP(qoiRoutine<P_V,P_M,Q_V,Q_M>, (void *) &valQoiRoutine_Data);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  cycle.valFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

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

  if (cycle.calFP().computeSolutionFlag() &&
      cycle.valFP().computeSolutionFlag()) {
    Q_V* epsilonVec = cycle.calFP().qoiRv().imageSet().vectorSpace().newVector(0.02);
    Q_V cdfDistancesVec(cycle.calFP().qoiRv().imageSet().vectorSpace().zeroVector());
    horizontalDistances(cycle.calFP().qoiRv().cdf(),
                        cycle.valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(cycle.valFP().qoiRv().cdf(),
                        cycle.calFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (cycle.env().rank() == 0) {
      std::cout << "For epsilonVec = "                             << *epsilonVec
                << ", cdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
  }

  // done
  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqBurgers_No_Model_Uncertainty' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_BURGERS_NO_MODEL_UNCERTAINTY_H__
