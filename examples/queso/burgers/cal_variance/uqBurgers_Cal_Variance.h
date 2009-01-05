/* uq/examples/queso/burgers/cal_variance/uqBurgers_Cal_Variance.h
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

#ifndef __UQ_BURGERS_CAL_VARIANCE_H__
#define __UQ_BURGERS_CAL_VARIANCE_H__

#define PI 3.14159265358979323846

// queso
#include <uqDefines.h>
#include <uqEnvironment.h>
#include <uqVectorSpace.h>
#include <uqVectorRV.h>
#include <uqAsciiTable.h>
#include <uqValidationCycle.h>
#include <uqGslVector.h>
#include <uqGslMatrix.h>

// system
#include <sys/time.h>

// gsl
#include <gsl/gsl_vector.h>

// burgers
#include <quesoInterfaceEMI.h>

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
			      const int N,
                              const char* fileName,
			      burgersQuesoInterface& interface);
 ~likelihoodRoutine_DataClass();

  // burgers/queso interface class
  burgersQuesoInterface& burgersInterface;

  // Experimental data
  const int nScenario; // number of data points

  // Scenario parameters
  int *nXloc; // number of x locations (1 for each scenario)
  double **xx; // x locations of data points
  double *Re; // Reynolds numbers of data points
  double *Freq; // Pressure grads of data points

  // Data values
  double **data; // either state (for calibration) or Reynolds stress (for validation)

};

// likelihoodRoutine_DataClass constructor: Read data, allocate/initialize memory for Burgers' solver
template<class P_V, class P_M>
likelihoodRoutine_DataClass<P_V,P_M>::likelihoodRoutine_DataClass(const uqBaseEnvironmentClass& env,
								  const int N,
								  const char* fileName,
								  burgersQuesoInterface& interface)
  :
  burgersInterface(interface),
  nScenario(N)
{
  std::cout << "Calling likelihoodRoutine_DataClass constructor..." << std::flush;

  // open file
  FILE *fp = fopen(fileName, "r");
  UQ_FATAL_TEST_MACRO((!fp), env.rank(),
		      "uqAppl(), in uqBurgers_EM_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		      "could not open data file");
  
  // get number of data points
  int nTmp;
  UQ_FATAL_TEST_MACRO( (fscanf(fp, "%d", &nTmp)!=1), env.rank(),
		       "uqAppl(), in uqBurgers_EM_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		       "could not read number of data points");

  UQ_FATAL_TEST_MACRO( (nTmp!=nScenario), env.rank(),
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		       "incorrect number of points in data file");
    
  nXloc = new int[nScenario];
  Re = new double[nScenario];
  Freq = new double[nScenario];

  xx = new double*[nScenario];
  data = new double*[nScenario];
  
  double tmpRe, tmpFreq, tmpXX, tmpU;

  // loop through scenarios
  for( int ii=0; ii<nScenario; ii++ ){
    fscanf(fp, "%d %lf %lf", &nTmp, &tmpRe, &tmpFreq);
    nXloc[ii] = nTmp;
    Re[ii] = tmpRe;
    Freq[ii] = tmpFreq;
    xx[ii] = new double[nTmp];
    data[ii] = new double[nTmp];

    for( int jj=0; jj<nTmp; jj++ ){
      fscanf(fp, "%lf %lf", &tmpXX, &tmpU);
      xx[ii][jj] = tmpXX;
      data[ii][jj] = tmpU;
    }
  }

  // close file
  fclose(fp);

  std::cout << "success.\n" << std::flush;
} // end constructor

// likelihoodRoutine_DataClass destructor: Frees memory allocated by constructor
template<class P_V, class P_M>
likelihoodRoutine_DataClass<P_V,P_M>::~likelihoodRoutine_DataClass()
{
  std::cout << "Calling likelihoodRoutine_DataClass destructor..." << std::flush;
  delete[] xx;
  delete[] Re;
  delete[] Freq;
  delete[] data;
  std::cout << "success.\n" << std::flush;
} // end destructor

// prototype for covariance, defined in .C file
int
computeLikelihood(likelihoodRoutine_DataClass<uqGslVectorClass, uqGslMatrixClass>* functionDataPtr, 
		  double modelVar, double expVar, double* misfit, double& nTwoLogLikelihood);

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
likelihoodRoutine(const P_V& paramValues, const void* functionDataPtr, P_V* gV, P_M* hM, P_V* hE)
{
  int ierr;

  // Experimental data
  const int nScenario = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->nScenario;
  const int *nXloc = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->nXloc;
  double **xx = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->xx;
  const double *Re = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->Re;
  const double *Freq = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->Freq;
  double **data = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->data;

  burgersQuesoInterface& burgersInterface = ((likelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->burgersInterface;

  // Parameters
  double kappa = paramValues[0];
  double modelVar = exp(paramValues[1])*exp(paramValues[1]); // the log of the standard dev is stored
  double myinadequacy[1] = {0.0}; // this is just a place holder... it ensures that the embedded model inadequacy is zero

  // Figure out total number of data points
  int nTot = 0;
  for( int ii=0; ii<nScenario; ii++ ) nTot += nXloc[ii];

  // Compute misfit vector
  double misfit[nTot];

  int index = 0;
  for( int ii=0; ii<nScenario; ii++ ){

    // Model results
    double model[nXloc[ii]];
    ierr = burgersInterface.solveForStateAtDataLocations(Re[ii], Freq[ii], kappa, myinadequacy,
							 nXloc[ii], xx[ii], model);

    if( ierr != 0 ){
      std::cerr << "WARNING: Solution did not converge... assuming params are non-physical.\n";
      // if solver was not successful, assume model is non-physical ==> Set likelihood very small
      return 1e30; // i.e. set -2*log(likelihood) very large
    }

    // misfit with data
    for( int jj=0; jj<nXloc[ii]; jj++ ){
      misfit[index + jj] = (data[ii][jj] - model[jj]);
    }
    index += nXloc[ii];
  }

  // Compute covariance matrix
  double expVar = 1e-10;
  double nTwoLogLikelihood=0;
  ierr = computeLikelihood((likelihoodRoutine_DataClass<P_V,P_M> *)functionDataPtr, modelVar, expVar, misfit,
			   nTwoLogLikelihood);
  UQ_FATAL_TEST_MACRO( (ierr!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_Cal_Variance (likelihoodRoutine)",
		       "failed during likelihood calculation" );

  gV = (P_V*) NULL;
  hM = (P_M*) NULL;
  hE = (P_V*) NULL;

  return nTwoLogLikelihood;
}


//********************************************************
// QoI function object 
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************

template<class P_V, class P_M, class Q_V, class Q_M>
struct
qoiRoutine_DataClass
{
  qoiRoutine_DataClass( uqVectorSpaceClass<Q_V,Q_M>& myQoiSpace,
			const uqBaseVectorRVClass<P_V,P_M>& myPostRv,
			burgersQuesoInterface& myBurgersInterface);
  ~qoiRoutine_DataClass();

  void qoiMonteCarloIntegrate();

  uqVectorSpaceClass<Q_V,Q_M>& qoiSpace;
  uqGenericVectorRVClass<Q_V,Q_M>* qoiRv;

  // Chain to use for quantity of interest
  const uqBaseVectorRVClass<P_V,P_M>& postRv;

  // Chains computed from posterior RV chain
  uqBaseVectorSequenceClass<Q_V, Q_M>* qoi_chain; // chain of qoi, evaled from chain of model params
  uqBaseVectorSequenceClass<Q_V, Q_M>* var_chain; // chain of qoi variance, evaled chain of model vars

  // Grid and values for CDF
  Q_V* numIntervalsVec;
  uqArrayOfOneDGridsClass<Q_V, Q_M>* cdfGrids;
  uqArrayOfOneDTablesClass<Q_V, Q_M>* cdfValues;

  // interfaces to burgers solver
  burgersQuesoInterface& burgersInterface;
};

// constructor
template<class P_V, class P_M, class Q_V, class Q_M>
qoiRoutine_DataClass<P_V, P_M, Q_V, Q_M>::qoiRoutine_DataClass( uqVectorSpaceClass<Q_V,Q_M>& myQoiSpace,
								const uqBaseVectorRVClass<P_V,P_M>& myPostRv,
								burgersQuesoInterface& myBurgersInterface)
  :
  qoiSpace(myQoiSpace),
  postRv(myPostRv),
  burgersInterface(myBurgersInterface)
{

  // call RV constructor (is this right??)
  qoiRv = new uqGenericVectorRVClass<Q_V, Q_M>("qoi_", this->qoiSpace);

  // allocate chains of length 0
  qoi_chain = new uqSequenceOfVectorsClass<Q_V, Q_M>(qoiRv->imageSet().vectorSpace(), 0, "qoi_chain"); 
  var_chain = new uqSequenceOfVectorsClass<Q_V, Q_M>(qoiRv->imageSet().vectorSpace(), 0, "var_chain");

  // cdf
  numIntervalsVec = qoiRv->imageSet().vectorSpace().newVector(250);
  
  cdfGrids  = new uqArrayOfOneDGridsClass<Q_V, Q_M>("qoi_cdf", qoiRv->imageSet().vectorSpace());

  Q_V *minQoi = qoiRv->imageSet().vectorSpace().newVector(5.0e-2);
  Q_V *maxQoi = qoiRv->imageSet().vectorSpace().newVector(7.5e-2);

  cdfGrids->setUniformGrids(*numIntervalsVec, *minQoi, *maxQoi);

  cdfValues = new uqArrayOfOneDTablesClass<Q_V, Q_M>("qoi_cdf", qoiRv->imageSet().vectorSpace());

} // end constructor

// destructor
template<class P_V, class P_M, class Q_V, class Q_M>
qoiRoutine_DataClass<P_V, P_M, Q_V, Q_M>::~qoiRoutine_DataClass()
{
  delete qoi_chain;
  delete var_chain;
  delete numIntervalsVec;
  delete cdfGrids;
  delete cdfValues;
} // end destructor


// Routine to compute CDF using Monte Carlo integration
// given properly filled qoi and mean chains
template<class P_V, class P_M, class Q_V, class Q_M>
void 
qoiRoutine_DataClass<P_V, P_M, Q_V, Q_M>::qoiMonteCarloIntegrate()
{
  int ierr;
  unsigned int seqSize = this->postRv.realizer().period();
  double zz, qoi, var, tmp, g;
  double myInadequacy = 0.0;
  //Q_V tmpQ(this->qoiRv->imageSet().vectorSpace().zeroVector());
  //Q_V tmpMean(this->qoiRv->imageSet().vectorSpace().zeroVector());

  const uqBaseOneDGridClass<double>& myGrid = this->cdfGrids->grid(0);
  std::vector<double> integral(myGrid.size(), 0.0);

  // From posterior chain, compute base qoi (parameter uncertainty only)
  // and qoi variance from model uncertainty
  printf("Computing qoi chain of length %d...\n", seqSize); fflush(stdout);

  // Compute output realizer: Monte Carlo approach
  this->qoi_chain->resizeSequence(seqSize);
  this->var_chain->resizeSequence(seqSize);

  P_V tmpV(this->postRv.imageSet().vectorSpace().zeroVector());
  Q_V tmpQ(this->qoiRv->imageSet().vectorSpace().zeroVector());
  Q_V tmpVar(this->qoiRv->imageSet().vectorSpace().zeroVector());

  for (unsigned int ii = 0; ii < seqSize; ii++) {
    this->postRv.realizer().realization(tmpV);

    // Compute model only qoi
    ierr = burgersInterface.solveForViscFluxAtOne(200.0, 1.0, tmpV[0], &myInadequacy, &(tmpQ[0]) );
    UQ_FATAL_TEST_MACRO( (ierr!=0), UQ_UNAVAILABLE_RANK,
			 "Burgers solver unsuccessful in qoiMonteCarloIntegrate",
			 "quitting" );

    this->qoi_chain->setPositionValues(ii,tmpQ);

    tmpVar[0] = exp(tmpV[1])*exp(tmpV[1])/(0.1*0.1); // 0.1 is length scale
    tmpVar[0] *= (5e-3)*(5e-3); // 5e-3 is viscosity at scenario of interest
    this->var_chain->setPositionValues(ii, tmpVar);

  }

  for( unsigned int ii=0; ii<myGrid.size(); ii++ ){  // loop through points where we want the cdf
    
    zz = myGrid[ii];
    integral[ii] = 0.0;

    for( unsigned int jj=0; jj<seqSize; jj++ ){ // loop through sequence

      // Compute qoi and qoi variance at this point on chain
      this->qoi_chain->getPositionValues(jj, tmpQ);
      qoi = tmpQ[0];
      
      this->var_chain->getPositionValues(jj, tmpVar);
      var = tmpVar[0];

      //std::cout << "var = " << std::scientific << std::setprecision(15) << var << "\n" << std::flush;
      
      // CORRECT
      // evaluate integrand
      tmp = (zz - qoi)/(sqrt(var*2.0));
      g = 0.5*(1.0 + erf(tmp));

//       // PARAMETER ONLY (i.e. var = 0)
//       if( zz > qoi ){
// 	g = 1.0;
//       } else{
// 	g = 0.0;
//       }

      // add to sum
      integral[ii] += g;
    }
    
    // set cdf value
    integral[ii] /= ((double)seqSize);
  }
  this->cdfValues->setOneDTable(0, integral);

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
    std::cout << "Beginning run of 'uqBurgers_EM_Model_Uncertainty' example\n"
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
					2, // # of rows FIX ME: (model parameter + model uncertainty variance)
					5,    // # of cols after 'parameter name': min + max + initial value for Markov chain
					NULL, // All extra columns are of 'double' type
					"param.tab");


  const EpetraExt::DistArray<std::string>& paramNames = paramTable.stringColumn(0);
  P_V paramMinValues    (paramTable.doubleColumn(1));
  P_V paramMaxValues    (paramTable.doubleColumn(2));
  P_V calInitialValues  (paramTable.doubleColumn(3));
  P_V paramExpValues    (paramTable.doubleColumn(4));
  P_V paramVarValues    (paramTable.doubleColumn(5));

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
                                      0,    // # of cols after 'parameter name': FIX ME: pass in qoi scenario???
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

//   // Set up inverse problem (i.e. prior and likelihood) FIX ME: Change to Gaussian prior RV??
//   uqUniformVectorRVClass<P_V,P_M> calPriorRv("cal_prior_", // Extra prefix before the default "rv_" prefix
//                                              paramDomain);

  // Set up inverse problem (i.e. prior and likelihood) FIX ME: Change to Gaussian prior RV??
  uqGaussianVectorRVClass<P_V,P_M> calPriorRv("cal_prior_", // Extra prefix before the default "rv_" prefix
					      paramDomain,
					      paramExpValues,
					      paramVarValues);

  // FIX ME: THIS IS THE WRONG BURGERS INTERFACE!!
  burgersQuesoInterface cal_burgersInterface(20, 1, 1, 1);// 20 solution modes, 1 inadequacy mode (not used here)
  //burgersQuesoInterface cal_burgersInterface(100, 1, 1, 1);// 20 solution modes, 1 inadequacy mode (not used here)

  // FIX ME: PASSING IN THE WRONG INTERFACE
  likelihoodRoutine_DataClass<P_V,P_M> calLikelihoodRoutine_Data(env, 2, "calibration_Re10-20_32pts.dat", 
								 cal_burgersInterface);
//   likelihoodRoutine_DataClass<P_V,P_M> calLikelihoodRoutine_Data(env, 3, "calibration_Re10-20-100_32pts.dat", 
// 								 cal_burgersInterface);


//   uqGenericScalarFunctionClass<P_V,P_M> calLikelihoodFunctionObj("cal_like_",
//                                                                  paramDomain,
//                                                                  likelihoodRoutine<P_V,P_M>,
//                                                                  NULL,
//                                                                  NULL,
//                                                                  (void *) &calLikelihoodRoutine_Data,
//                                                                  true);

  uqGenericScalarFunctionClass<P_V,P_M> calLikelihoodFunctionObj("cal_like_",
                                                                 paramDomain,
                                                                 likelihoodRoutine<P_V,P_M>,
                                                                 (void *) &calLikelihoodRoutine_Data,
                                                                 true);

  // Set up calibration inverse problem
  cycle.setCalIP(calPriorRv, calLikelihoodFunctionObj);


  // Solve calibration inverse problem = set 'pdf' and 'realizer' of 'postRv' (using Markov Chain)
  P_M* calProposalCovMatrix = 
    cycle.calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(cycle.calIP().priorRv().pdf().domainVarVector(),
								      calInitialValues);

  cycle.calIP().solveWithBayesMarkovChain(calInitialValues, calProposalCovMatrix);
  delete calProposalCovMatrix;

  
  // Forward propagation problem
  burgersQuesoInterface prop_burgersInterface(200, 1, 1, 1);//
  qoiRoutine_DataClass<P_V, P_M, Q_V, Q_M> calQoiRoutine_Data(qoiSpace, cycle.calIP().postRv(), prop_burgersInterface);
  calQoiRoutine_Data.qoiMonteCarloIntegrate();

  uqBaseVectorCdfClass<Q_V,Q_M>* calQoiCdf = new uqSampledVectorCdfClass<Q_V,Q_M>("cal_fp_qoi_",
										  calQoiRoutine_Data.cdfGrids[0],
										  calQoiRoutine_Data.cdfValues[0]);

  calQoiRoutine_Data.qoiRv->setCdf(*calQoiCdf);

  // Write output file
  if (env.rank() == 0) {
    std::cout << "Opening output file '" << "calOutput.m"
	      << "' for propagation problem with problem with prefix = " << "cal_fp"
	      << std::endl;
  }
  
  // Open file
  std::ofstream* ofs = new std::ofstream("calOutput.m", std::ofstream::out | std::ofstream::in | std::ofstream::ate);
  if ((ofs            == NULL ) ||
      (ofs->is_open() == false)) {
    delete ofs;
    ofs = new std::ofstream("calOutput.m", std::ofstream::out | std::ofstream::trunc);
  }
  UQ_FATAL_TEST_MACRO((ofs && ofs->is_open()) == false,
		      env.rank(),
		      "uqPropagProblem<P_V,P_M,Q_V,Q_M>::solveWithBayesMarkovChain()",
		      "failed to open file");
  
  *ofs << calQoiRoutine_Data.qoiRv->cdf();
  
  // Close file
  ofs->close();
  delete ofs;
  if (env.rank() == 0) {
    std::cout << "Closed output file '" << "calOutput.m"
	      << "' for propagation problem with problem with prefix = " << "cal_fp"
	      << std::endl;
  }


  // done with calibration phase
  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending 'calibration stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'calibration stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }


//   //******************************************************
//   // Burgers' example: validation phase
//   //******************************************************
//   iRC = gettimeofday(&timevalRef, NULL);
//   if (env.rank() == 0) {
//     std::cout << "Beginning 'validation stage' at " << ctime(&timevalRef.tv_sec)
//               << std::endl;
//   }

//   burgersQuesoInterface val_burgersInterface(100, 1, 1, 1);// 20 solution modes, 1 inadequacy mode (not used here)

//   likelihoodRoutine_DataClass<P_V,P_M> valLikelihoodRoutine_Data(env, 1, "validation_Re100_32pts.dat", 
// 								 val_burgersInterface);


//   uqGenericScalarFunctionClass<P_V,P_M> valLikelihoodFunctionObj("val_like_",
//                                                                  paramDomain,
//                                                                  likelihoodRoutine<P_V,P_M>,
//                                                                  NULL,
//                                                                  NULL,
//                                                                  (void *) &valLikelihoodRoutine_Data,
//                                                                  true);
//   // Set up validation inverse problem
//   cycle.setValIP(valLikelihoodFunctionObj);


//   // Solve valibration inverse problem = set 'pdf' and 'realizer' of 'postRv' (using Markov Chain)
//   P_M* valProposalCovMatrix = 
//     cycle.calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(cycle.calIP().postRv().realizer().imageVarVector(),
// 								      cycle.calIP().postRv().realizer().imageExpVector() );

//   cycle.valIP().solveWithBayesMarkovChain(cycle.calIP().postRv().realizer().imageExpVector(),
//                                           valProposalCovMatrix);
//   delete valProposalCovMatrix;

  
//   // Forward propagation problem
//   //burgersQuesoInterface prop_burgersInterface(200, 1, 1, 1);// already instantiated
//   qoiRoutine_DataClass<P_V, P_M, Q_V, Q_M> valQoiRoutine_Data(qoiSpace, cycle.valIP().postRv(), prop_burgersInterface);
//   valQoiRoutine_Data.qoiMonteCarloIntegrate();

//   uqBaseVectorCdfClass<Q_V,Q_M>* valQoiCdf = new uqSampledVectorCdfClass<Q_V,Q_M>("val_fp_qoi_",
// 										  valQoiRoutine_Data.cdfGrids[0],
// 										  valQoiRoutine_Data.cdfValues[0]);

//   valQoiRoutine_Data.qoiRv->setCdf(*valQoiCdf);

//   // Write output file
//   if (env.rank() == 0) {
//     std::cout << "Opening output file '" << "valOutput.m"
// 	      << "' for propagation problem with problem with prefix = " << "val_fp"
// 	      << std::endl;
//   }
  
//   // Open file
//   ofs = new std::ofstream("valOutput.m", std::ofstream::out | std::ofstream::in | std::ofstream::ate);
//   if ((ofs            == NULL ) ||
//       (ofs->is_open() == false)) {
//     delete ofs;
//     ofs = new std::ofstream("valOutput.m", std::ofstream::out | std::ofstream::trunc);
//   }
//   UQ_FATAL_TEST_MACRO((ofs && ofs->is_open()) == false,
// 		      env.rank(),
// 		      "uqPropagProblem<P_V,P_M,Q_V,Q_M>::solveWithBayesMarkovChain()",
// 		      "failed to open file");
  
//   *ofs << valQoiRoutine_Data.qoiRv->cdf();
  
//   // Close file
//   ofs->close();
//   delete ofs;
//   if (env.rank() == 0) {
//     std::cout << "Closed output file '" << "valOutput.m"
// 	      << "' for propagation problem with problem with prefix = " << "val_fp"
// 	      << std::endl;
//   }


//   // done with valibration phase
//   iRC = gettimeofday(&timevalNow, NULL);
//   if (env.rank() == 0) {
//     std::cout << "Ending 'valibration stage' at " << ctime(&timevalNow.tv_sec)
//               << "Total 'valibration stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
//               << " seconds"
//               << std::endl;
//   }


//   // FIX ME: Add comparison phase



  return;
}
#endif // __UQ_BURGERS_CAL_VARIANCE_H__
