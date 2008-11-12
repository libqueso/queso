/* uq/examples/queso/burgers/gp_model_uncertainty/uqBurgers_GP_Model_Uncertainty.h
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


// A struct used by the likelihood data
// This struct contains experimental data and
// a covariance matrix specifying model uncertainty
// covariance at the data points
struct
inverseProblem_DataClass
{
  inverseProblem_DataClass(const uqEnvironmentClass& env,
			   double Re,
			   const char* fileName);
  ~inverseProblem_DataClass();
  
  // Experimental data
  int nDataPoints;       // number of data points
  double Re;             // Reynolds number for data
  double *dataLocations; // x locations where state was measured
  double *dataValues;    // measured values of averaged state at dataLocations
  
  // Gaussian process model uncertainty
  double *cholK; // GP covariance evaluated at data points
};
  
// inverseProblem_DataClass constructor
inverseProblem_DataClass::inverseProblem_DataClass(const uqEnvironmentClass& env,
						   double Re,
						   const char* fileName)
  :
  nDataPoints(0),
  dataLocations(0),
  dataValues(0),
  cholK(0)
{
  // Set Re
  this->Re = Re;

  // open data file
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

  // done with data file
  fclose(fp);

  // allocate storage for covariance and initialize to identity
  cholK = new double[nDataPoints*nDataPoints];
  for( int ii=0; ii<nDataPoints; ii++ ){
    for( int jj=0; jj<nDataPoints; jj++ ){
      cholK[nDataPoints*ii+jj] = 0.0;
    }
    cholK[nDataPoints*ii+ii] = 1.0;
  }

} // end constructor

// inverseProblem_DataClass destructor
inverseProblem_DataClass::~inverseProblem_DataClass()
{
  delete[] dataLocations;
  delete[] dataValues;
  delete[] cholK;
} // end destructor
  

// Declaration of function to compute GP covariance matrix for calibration phase
// (defined in uqBurgers_GP_Model_Uncertainty.C)
int
calibrationCovarianceChol(const double sig2, const double ellx, inverseProblem_DataClass *cal_data);

// Declaration of function to compute GP covariance matrix for validation phase
// (defined in uqBurgers_GP_Model_Uncertainty.C)
int
validationCovarianceChol(const double sig2, const double ellx, const double ellRe,
			 const inverseProblem_DataClass *cal_data, inverseProblem_DataClass *val_data);


//********************************************************
// Likelihood function object for both inverse problems of the validation cycle.
// A likelihood function object is provided by user and is called by the UQ library.
// This likelihood function object consists of data and routine.
//********************************************************

// The (user defined) data class for the data needed by the 
// (user defined) likelihood routine
struct 
likelihoodRoutine_DataClass
{
  likelihoodRoutine_DataClass(const uqEnvironmentClass& env,
			      double cal_Re,
			      double val_Re,
                              const char* calFileName,
			      const char* valFileName);
 ~likelihoodRoutine_DataClass();

  // Experimental data and model uncertainty
  bool calPhase;                             // if true, in calibration phase
  inverseProblem_DataClass *calibrationData; // calibration data and covariance
  inverseProblem_DataClass *validationData;  // validation data and covariance

  // For Burgers' solver
  quadBasis *pQB; // quad points, weights, and basis evaluation
  gsl_vector *U;  // solution (re-used after each solve at IC for next solve to hopefully decrease Newton iterations)
};

// Declaration of function to compute likelihood from misfit
// (defined in uqBurgers_GP_Model_Uncertainty_gsl.C)
void 
negTwoLogLikelihood(const int N, const double *cholK, const double *misfit, 
		    double *result);

// Declaration of function to compute GP mean for validation phase (i.e. posterior mean from calibration)
// (defined in uqBurgers_GP_Model_Uncertainty_gsl.C)
int
validationMean(const double sig2, const double ellx, const double ellRe, const double kappa,
	       likelihoodRoutine_DataClass *data, double *mean);


// likelihoodRoutine_DataClass constructor: Read data, allocate/initialize memory for Burgers' solver
likelihoodRoutine_DataClass::likelihoodRoutine_DataClass(const uqEnvironmentClass &env, 
							 double cal_Re, 
							 double val_Re,
							 const char *calFileName,
							 const char *valFileName)
  :
  pQB(0),
  U(0)
{
  // allocation and initialize calibration and/or validation data
  calPhase = false;
  calibrationData = (inverseProblem_DataClass *)NULL;
  validationData  = (inverseProblem_DataClass *)NULL;

  if( calFileName ){
    calPhase = true;
    calibrationData = new inverseProblem_DataClass(env, cal_Re, calFileName);
  }

  if( valFileName ){
    calPhase = false;  // ONLY have validation data for validation phase
    validationData = new inverseProblem_DataClass(env, val_Re, valFileName);
  }

  // allocate and initialize quadrature points
  UQ_FATAL_TEST_MACRO( (evaluateQuadratureAndBasisForResidual(100, &pQB)!=0), env.rank(),
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		       "failed while allocating/computing quadrature points" );
		       
  // allocate and initialize solution (to zero) 
  U = gsl_vector_calloc(100);

} // end constructor


// likelihoodRoutine_DataClass destructor: Frees memory allocated by constructor
likelihoodRoutine_DataClass::~likelihoodRoutine_DataClass()
{
  delete calibrationData;
  delete validationData;
  freeQuadratureAndBasisForResidual(pQB);
  gsl_vector_free(U);
} // end destructor


// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
likelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  // FIXME: This routine is ugly

  int nDataPoints;
  double Re, *xx, *ue, *cholK;

  int cal_nDataPoints;
  double cal_Re, *cal_xx, *cal_ue, *cal_cholK;

  int val_nDataPoints;
  double val_Re, *val_xx, *val_ue, *val_cholK;

  bool calPhase = ((likelihoodRoutine_DataClass *) functionDataPtr)->calPhase;

  // Experimental data
  cal_nDataPoints = ((likelihoodRoutine_DataClass *) functionDataPtr)->calibrationData->nDataPoints;
  cal_Re = ((likelihoodRoutine_DataClass *) functionDataPtr)->calibrationData->Re;
  cal_xx = ((likelihoodRoutine_DataClass *) functionDataPtr)->calibrationData->dataLocations;
  cal_ue = ((likelihoodRoutine_DataClass *) functionDataPtr)->calibrationData->dataValues;

  cal_cholK = ((likelihoodRoutine_DataClass *) functionDataPtr)->calibrationData->cholK;
    
  // Gaussian process model uncertainty
  if( !calPhase ){
    // Experimental data
    val_nDataPoints = ((likelihoodRoutine_DataClass *) functionDataPtr)->validationData->nDataPoints;
    val_Re = ((likelihoodRoutine_DataClass *) functionDataPtr)->validationData->Re;
    val_xx = ((likelihoodRoutine_DataClass *) functionDataPtr)->validationData->dataLocations;
    val_ue = ((likelihoodRoutine_DataClass *) functionDataPtr)->validationData->dataValues;
    
    // Gaussian process model uncertainty
    val_cholK = ((likelihoodRoutine_DataClass *) functionDataPtr)->validationData->cholK;
  }

  // Solver data
  quadBasis *pQB = ((likelihoodRoutine_DataClass *) functionDataPtr)->pQB;
  gsl_vector *U = ((likelihoodRoutine_DataClass *) functionDataPtr)->U;

  // Model parameter
  double kappa = paramValues[0];

  // choose appropriate data
  if( calPhase ){
    nDataPoints = cal_nDataPoints;
    Re = cal_Re;
    xx = cal_xx;
    ue = cal_ue;
    cholK = cal_cholK;
  } else{
    nDataPoints = val_nDataPoints;
    Re = val_Re;
    xx = val_xx;
    ue = val_ue;
    cholK = val_cholK;
  }

  // Model data
  double um[nDataPoints];

  // compute GP mean
  double mean[nDataPoints];
  for( int ii=0; ii<nDataPoints; ii++ ) mean[ii] = 0.0; // for calibration phase, model inadequacy mean is zero
  
  if( !calPhase ){ // for validation phase, model inadequacy mean is posterior from calibration
    int ierr = validationMean(1e-4, 0.1, 0.5, kappa, (likelihoodRoutine_DataClass *) functionDataPtr, mean);
    UQ_FATAL_TEST_MACRO( (ierr!=0), UQ_UNAVAILABLE_RANK,
			 "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine)",
			 "failed while computing validation phase GP mean" );
  }

  // Evaluate model
  int ierr = solveForStateAtXLocations(xx, nDataPoints, 1.0/Re, kappa, U, pQB, um);
  UQ_FATAL_TEST_MACRO( (ierr!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine)",
		       "failed solving Burgers' eqn" );

  // Compute misfit vector
  double misfit[nDataPoints];
  for( int ii=0; ii<nDataPoints; ii++ ) misfit[ii] = (ue[ii] - um[ii] - mean[ii]);

  // Compute -2*log(likelihood)
  double resultValue;
  negTwoLogLikelihood(nDataPoints, cholK, misfit, &resultValue);

  return resultValue;
}


//********************************************************
// QoI function object 
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data type for the data needed by the (user defined) qoi routine
template<class Q_V, class Q_M>
struct
qoiRoutine_DataClass
{
  qoiRoutine_DataClass( uqVectorSpaceClass<Q_V,Q_M>& myQoiSpace );
  ~qoiRoutine_DataClass();

  void qoiMonteCarloIntegrate();

  uqVectorSpaceClass<Q_V,Q_M>& qoiSpace;
  uqGenericVectorRVClass<Q_V,Q_M>* qoiRv;

  // Data for computing qoi CDF
  uqBaseVectorSequenceClass<Q_V, Q_M>* qoi_chain; // chain of qoi, evaled from chain of params
  uqBaseVectorSequenceClass<Q_V, Q_M>* mean_chain; // chain of GP posterior mean for qoi, evaled from chain of params
  double *variance; // GP posterior variance for qoi (independent of params)

  // Grid and values for CDF
  Q_V* numIntervalsVec;
  uqArrayOfOneDGridsClass<Q_V, Q_M>* cdfGrids;
  uqArrayOfOneDTablesClass<Q_V, Q_M>* cdfValues;

  // For the Burgers solver
  quadBasis *pQB; // quad points, weights, and basis evaluation
};

// constructor
template<class Q_V, class Q_M>
qoiRoutine_DataClass<Q_V, Q_M>::qoiRoutine_DataClass( uqVectorSpaceClass<Q_V,Q_M>& myQoiSpace )
  :
  qoiSpace(myQoiSpace)
{

  // call RV constructor (is this right??)
  qoiRv = new uqGenericVectorRVClass<Q_V, Q_M>("qoi_", this->qoiSpace);

  // allocate chains of length 0
  qoi_chain  = new uqSequenceOfVectorsClass<Q_V, Q_M>(qoiRv->imageSpace(), 0, "qoi_chain"); 
  mean_chain = new uqSequenceOfVectorsClass<Q_V, Q_M>(qoiRv->imageSpace(), 0, "mean_chain");

  variance = new double[1];

  numIntervalsVec = qoiRv->imageSpace().newVector(250);
  
  cdfGrids  = new uqArrayOfOneDGridsClass<Q_V, Q_M>("qoi_cdf", qoiRv->imageSpace());

  Q_V *minQoi = qoiRv->imageSpace().newVector(5.0e-2);
  Q_V *maxQoi = qoiRv->imageSpace().newVector(7.5e-2);

  cdfGrids->setUniformGrids(*numIntervalsVec, *minQoi, *maxQoi);

  cdfValues = new uqArrayOfOneDTablesClass<Q_V, Q_M>("qoi_cdf", qoiRv->imageSpace());

  // allocate and initialize quadrature points
  UQ_FATAL_TEST_MACRO( (evaluateQuadratureAndBasisForResidual(100, &pQB)!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine_DataClass constructor)",
		       "failed while allocating/computing quadrature points" );
} // end constructor

// destructor
template<class Q_V, class Q_M>
qoiRoutine_DataClass<Q_V, Q_M>::~qoiRoutine_DataClass()
{
  delete qoi_chain;
  delete mean_chain;
  delete[] variance;
  delete numIntervalsVec;
  delete cdfGrids;
  delete cdfValues;
  freeQuadratureAndBasisForResidual(pQB);
} // end destructor


// Routine to compute CDF using Monte Carlo integration
// given properly filled qoi and mean chains
template<class Q_V, class Q_M>
void 
qoiRoutine_DataClass<Q_V, Q_M>::qoiMonteCarloIntegrate()
{
  unsigned int seqSize = this->qoi_chain->sequenceSize();
  double zz, qoi, mean, tmp, g;
  Q_V tmpQ(this->qoiRv->imageSpace().zeroVector());
  Q_V tmpMean(this->qoiRv->imageSpace().zeroVector());

  const uqBaseOneDGridClass<double>& myGrid = this->cdfGrids->grid(0);
  std::vector<double> integral(myGrid.size(), 0.0);

  for( unsigned int ii=0; ii<myGrid.size(); ii++ ){  // loop through points where we want the cdf
    
    zz = myGrid[ii];
    integral[ii] = 0.0;

    for( unsigned int jj=0; jj<seqSize; jj++ ){ // loop through sequence

      // extract qoi and mean error from sequences previously computed
      this->qoi_chain->getPositionValues(jj, tmpQ);
      qoi = tmpQ[0];
      
      mean_chain->getPositionValues(jj, tmpMean);
      mean = tmpMean[0];
      
      // evaluate integrand
      tmp = (zz - qoi - mean)/(sqrt(this->variance[0]*2.0));
      g = 0.5*(1.0 + erf(tmp));

      // add to sum
      integral[ii] += g;
    }
    
    // set cdf value
    integral[ii] /= ((double)seqSize);
  }
  this->cdfValues->setOneDTable(0, integral);

  return;
}


int
computeStage1PropagationMean(const int N1, const double *xLoc1, const double Re1,
			     const int N3, const double *xLoc3, const double Re3,
			     const double sig2, const double ellx, const double ellRe,
			     const double *s1_misfit, const double *cholK, double *mean);

int
computeStage1PropagationVariance(const int N1, const double *xLoc1, const double Re1,
				 const int N3, const double *xLoc3, const double Re3,
				 const double sig2, const double ellx, const double ellRe,
				 const double *cholK, double *variance);

// Routine to compute qoi and mean chains for calibration phase
template<class P_V,class P_M,class Q_V, class Q_M>
void
calibrationQoiMeanChains(const uqBaseVectorRVClass<P_V,P_M>& calPostRv,
			 likelihoodRoutine_DataClass& calLikelihoodRoutine_Data,
			 qoiRoutine_DataClass<Q_V, Q_M>& qoiRoutine_Data)
{
  int ierr;
  unsigned int seqSize = calPostRv.realizer().period();
  quadBasis *pQB = qoiRoutine_Data.pQB;

  int N1 = calLikelihoodRoutine_Data.calibrationData->nDataPoints;
  double Re = calLikelihoodRoutine_Data.calibrationData->Re;
  double *xLoc1 = calLikelihoodRoutine_Data.calibrationData->dataLocations;
  double *ue1 = calLikelihoodRoutine_Data.calibrationData->dataValues;
  double xLoc3 = 1.0;
  double *cholK = calLikelihoodRoutine_Data.calibrationData->cholK;
  double variance;

  gsl_vector *U1 = calLikelihoodRoutine_Data.U;
  double um1[N1];


  // Compute variance of conditional dist of d(eps)/dx (it is independent of kappa)
  ierr = computeStage1PropagationVariance(N1, xLoc1, Re, 1, &xLoc3, 200, 1e-4, 0.1, 0.5, cholK, &variance);
  //
  variance *= (5e-3*5e-3); // FIX ME... DO THE RIGHT THING INSIDE THE FUNCTION!!!
  qoiRoutine_Data.variance[0] = variance;

  printf("Computing qoi chain of length %d...\n", seqSize); fflush(stdout);

  // Compute output realizer: Monte Carlo approach
  qoiRoutine_Data.qoi_chain->resizeSequence(seqSize);
  qoiRoutine_Data.mean_chain->resizeSequence(seqSize);

  P_V tmpV(calPostRv.imageSpace().zeroVector());
  Q_V tmpQ(qoiRoutine_Data.qoiRv->imageSpace().zeroVector());
  Q_V tmpMean(qoiRoutine_Data.qoiRv->imageSpace().zeroVector());

  for (unsigned int i = 0; i < seqSize; ++i) {
    calPostRv.realizer().realization(tmpV);

    // compute quantity of interest chain
    ierr = computeGradientAtOne(tmpV[0], pQB, &(tmpQ[0]));
    if( ierr != 0 ){
      printf("WARNING: Burgers solver was not successful in propagQoiRoutine!\n");
      printf("         Results are probably meaningless.\n");
      fflush(stdout);
    }
    qoiRoutine_Data.qoi_chain->setPositionValues(i,tmpQ);

    // compute mean value of d(eps)/dx for this kappa
    // evaluate misfit for stage 1 scenario
    ierr = solveForStateAtXLocations(xLoc1, N1, 1.0/Re, tmpV[0], U1, pQB, um1);
    if( ierr != 0 ){
      printf("WARNING: Burgers solver was not successful!\n");
      printf("         Results are probably meaningless.\n");
      fflush(stdout);
    }

    double s1_misfit[N1];
    for( int jj=0; jj<N1; jj++ ) s1_misfit[jj] = (ue1[jj] - um1[jj]);

    ierr = computeStage1PropagationMean(N1, xLoc1, Re, 1, &xLoc3, 200, 1e-4, 0.1, 0.5, 
					s1_misfit, cholK, &(tmpMean[0]));

    tmpMean[0] *= -5e-3;
    qoiRoutine_Data.mean_chain->setPositionValues(i, tmpMean);

  }
  
  return;
}

int
computeStage2PropagationMean(const int N1, const double *xLoc1, const double Re1,
			     const int N2, const double *xLoc2, const double Re2,
			     const int N3, const double *xLoc3, const double Re3,
			     const double sig2, const double ellx, const double ellRe,
			     const double *cholK1, const double *cholK2, double *s2_misfit, double mean1,
			     double *mean2);

int
computeStage2PropagationVariance(const int N1, const double *xLoc1, const double Re1,
				 const int N2, const double *xLoc2, const double Re2,
				 const int N3, const double *xLoc3, const double Re3,
				 const double sig2, const double ellx, const double ellRe,
				 const double *cholK1, const double *cholK2, const double variance1,
				 double *variance2);

// Routine to compute qoi and mean chains for calibration phase
template<class P_V,class P_M,class Q_V, class Q_M>
void
validationQoiMeanChains(const uqBaseVectorRVClass<P_V,P_M>& valPostRv,
			likelihoodRoutine_DataClass& valLikelihoodRoutine_Data,
			qoiRoutine_DataClass<Q_V, Q_M>& qoiRoutine_Data)
{
  int ierr;
  unsigned int seqSize = valPostRv.realizer().period();
  quadBasis *pQB = qoiRoutine_Data.pQB;

  int N1 = valLikelihoodRoutine_Data.calibrationData->nDataPoints;
  double Re1 = valLikelihoodRoutine_Data.calibrationData->Re;
  double *xLoc1 = valLikelihoodRoutine_Data.calibrationData->dataLocations;
  double *ue1 = valLikelihoodRoutine_Data.calibrationData->dataValues;
  double *cholK1 = valLikelihoodRoutine_Data.calibrationData->cholK;

  int N2 = valLikelihoodRoutine_Data.validationData->nDataPoints;
  double Re2 = valLikelihoodRoutine_Data.validationData->Re;
  double *xLoc2 = valLikelihoodRoutine_Data.validationData->dataLocations;
  double *ue2 = valLikelihoodRoutine_Data.validationData->dataValues;
  double *cholK2 = valLikelihoodRoutine_Data.validationData->cholK;

  double xLoc3 = 1.0;
  double variance1, variance2;

  gsl_vector *U1 = valLikelihoodRoutine_Data.U;
  double um1[N1], um2[N2], cal_qoi_mean;


  // Compute variance of conditional dist of d(eps)/dx (it is independent of kappa)
  ierr = computeStage1PropagationVariance(N1, xLoc1, Re1, 1, &xLoc3, 200, 1e-4, 0.1, 0.5, cholK1, &variance1);

  ierr = computeStage2PropagationVariance(N1, xLoc1, Re1, N2, xLoc2, Re2, 1, &xLoc3, 200.0,
					  1e-4, 0.1, 0.5, cholK1, cholK2, variance1, &variance2);
  //
  variance2 *= (5e-3*5e-3); // FIX ME... DO THE RIGHT THING INSIDE THE FUNCTION!!!
  qoiRoutine_Data.variance[0] = variance2;

  printf("Computing qoi chain of length %d...\n", seqSize); fflush(stdout);

  // Compute output realizer: Monte Carlo approach
  qoiRoutine_Data.qoi_chain->resizeSequence(seqSize);
  qoiRoutine_Data.mean_chain->resizeSequence(seqSize);

  P_V tmpV(valPostRv.imageSpace().zeroVector());
  Q_V tmpQ(qoiRoutine_Data.qoiRv->imageSpace().zeroVector());
  Q_V tmpMean(qoiRoutine_Data.qoiRv->imageSpace().zeroVector());

  for (unsigned int i = 0; i < seqSize; ++i) {
    valPostRv.realizer().realization(tmpV);

    // compute quantity of interest chain
    ierr = computeGradientAtOne(tmpV[0], pQB, &(tmpQ[0]));
    if( ierr != 0 ){
      printf("WARNING: Burgers solver was not successful in propagQoiRoutine!\n");
      printf("         Results are probably meaningless.\n");
      fflush(stdout);
    }
    qoiRoutine_Data.qoi_chain->setPositionValues(i,tmpQ);

    // compute mean value of d(eps)/dx for this kappa
    // evaluate misfit for stage 1 scenario
    ierr = solveForStateAtXLocations(xLoc1, N1, 1.0/Re1, tmpV[0], U1, pQB, um1);
    if( ierr != 0 ){
      printf("WARNING: Burgers solver was not successful!\n");
      printf("         Results are probably meaningless.\n");
      fflush(stdout);
    }

    double s1_misfit[N1];
    for( int jj=0; jj<N1; jj++ ) s1_misfit[jj] = (ue1[jj] - um1[jj]);

    ierr = computeStage1PropagationMean(N1, xLoc1, Re1, 1, &xLoc3, 200, 1e-4, 0.1, 0.5, 
					s1_misfit, cholK1, &cal_qoi_mean);


    // compute stage 2 mean value of d(eps)/dx for this kappa
    ierr = solveForStateAtXLocations(xLoc2, N2, 1.0/Re2, tmpV[0], U1, pQB, um2);
    if( ierr != 0 ){
      printf("WARNING: Burgers solver was not successful!\n");
      printf("         Results are probably meaningless.\n");
      fflush(stdout);
    }

    double cal_mean[N2];
    ierr = validationMean(1e-4, 0.1, 0.5, tmpV[0], &valLikelihoodRoutine_Data, cal_mean);

    double s2_misfit[N2];
    for( int jj=0; jj<N2; jj++ ) s2_misfit[jj] = (ue2[jj] - um2[jj] - cal_mean[jj]);

    ierr = computeStage2PropagationMean(N1, xLoc1, Re1, N2, xLoc2, Re2, 1, &xLoc3, 200.0,
					1e-4, 0.1, 0.5, cholK1, cholK2, s2_misfit, cal_qoi_mean, &(tmpMean[0]));

    tmpMean[0] *= -5e-3;
    qoiRoutine_Data.mean_chain->setPositionValues(i, tmpMean);

  }
  
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

  likelihoodRoutine_DataClass calLikelihoodRoutine_Data(env, 10.0, 100.0, "calibration_data.dat", NULL);

  // Set GP covariance for calibration phase (from prior)
  ierr = calibrationCovarianceChol(1e-4, 0.1, calLikelihoodRoutine_Data.calibrationData);
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

  
  // Forward propagation problem
  qoiRoutine_DataClass<Q_V, Q_M> calQoiRoutine_Data(qoiSpace);
  calibrationQoiMeanChains(cycle.calIP().postRv(), calLikelihoodRoutine_Data, calQoiRoutine_Data);
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


  //******************************************************
  // Burgers' example: validation phase
  //******************************************************
  iRC = gettimeofday(&timevalRef, NULL);
  if (env.rank() == 0) {
    std::cout << "Beginning 'validation stage' at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Set up inverse problem
  likelihoodRoutine_DataClass valLikelihoodRoutine_Data(env, 10.0, 100.0, "calibration_data.dat", "validation_data.dat");

  // Set GP covariance for calibration phase (from prior) (b/c this is required for validation phase)
  ierr = calibrationCovarianceChol(1e-4, 0.1, valLikelihoodRoutine_Data.calibrationData);
  UQ_FATAL_TEST_MACRO( (ierr!=0), env.rank(),
		       "uqAppl(), in uqBurgers_GP_Model_Uncertainty",
		       "failed while computing calibration covariance matrix" );

  // Set GP covariance for validation phase (from calibration posterior)
  ierr = validationCovarianceChol(1e-4, 0.1, 0.5, valLikelihoodRoutine_Data.calibrationData, 
				  valLikelihoodRoutine_Data.validationData);
  UQ_FATAL_TEST_MACRO( (ierr!=0), env.rank(),
		       "uqAppl(), in uqBurgers_GP_Model_Uncertainty",
		       "failed while computing validation covariance matrix" );

  cycle.setValIP(likelihoodRoutine<P_V,P_M>, (void *) &valLikelihoodRoutine_Data,
                 true); // the likelihood routine computes [-2.*ln(Likelihood)]


  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  // Use 'realizer()' because the posterior rv was computed with Markov Chain
  // Use these calibration mean values as the initial values
  P_M* valProposalCovMatrix =
    cycle.calIP().postRv().imageSpace().newGaussianMatrix(cycle.calIP().postRv().realizer().imageVarianceValues(),
							  cycle.calIP().postRv().realizer().imageExpectedValues());

  cycle.valIP().solveWithBayesMarkovChain(cycle.calIP().postRv().realizer().imageExpectedValues(),
                                          *valProposalCovMatrix,
                                          NULL); // use default kernel from library
  delete valProposalCovMatrix;

  // Forward propagation problem
  qoiRoutine_DataClass<Q_V, Q_M> valQoiRoutine_Data(qoiSpace);
  validationQoiMeanChains(cycle.valIP().postRv(), valLikelihoodRoutine_Data, valQoiRoutine_Data);
  valQoiRoutine_Data.qoiMonteCarloIntegrate();

  uqBaseVectorCdfClass<Q_V,Q_M>* valQoiCdf = new uqSampledVectorCdfClass<Q_V,Q_M>("val_fp_qoi_",
										  valQoiRoutine_Data.cdfGrids[0],
										  valQoiRoutine_Data.cdfValues[0]);

  valQoiRoutine_Data.qoiRv->setCdf(*valQoiCdf);

  // Write output file
  if (env.rank() == 0) {
    std::cout << "Opening output file '" << "valOutput.m"
	      << "' for propagation problem with problem with prefix = " << "val_fp"
	      << std::endl;
  }
  
  // Open file
  ofs = new std::ofstream("valOutput.m", std::ofstream::out | std::ofstream::in | std::ofstream::ate);
  if ((ofs            == NULL ) ||
      (ofs->is_open() == false)) {
    delete ofs;
    ofs = new std::ofstream("valOutput.m", std::ofstream::out | std::ofstream::trunc);
  }
  UQ_FATAL_TEST_MACRO((ofs && ofs->is_open()) == false,
		      env.rank(),
		      "uqPropagProblem<P_V,P_M,Q_V,Q_M>::solveWithBayesMarkovChain()",
		      "failed to open file");
  
  *ofs << valQoiRoutine_Data.qoiRv->cdf();
  
  // Close file
  ofs->close();
  delete ofs;
  if (env.rank() == 0) {
    std::cout << "Closed output file '" << "valOutput.m"
	      << "' for propagation problem with problem with prefix = " << "val_fp"
	      << std::endl;
  }

  // done with validation phase
  iRC = gettimeofday(&timevalNow, NULL);
  if (env.rank() == 0) {
    std::cout << "Ending 'validation stage' at " << ctime(&timevalNow.tv_sec)
              << "Total 'validation stage' run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }


  //******************************************************
  // Burgers' example: comparison phase
  //******************************************************

  // FIXME: Add comparison
  Q_V* epsilonVec = calQoiRoutine_Data.qoiRv->imageSpace().newVector(0.05);
  Q_V cdfDistancesVec(calQoiRoutine_Data.qoiRv->imageSpace().zeroVector());
  horizontalDistances(calQoiRoutine_Data.qoiRv->cdf(),
		      valQoiRoutine_Data.qoiRv->cdf(),
		      *epsilonVec,
		      cdfDistancesVec);
  if (env.rank() == 0) {
    std::cout << "For epsilonVec = "    << *epsilonVec
	      << ", cdfDistancesVec = " << cdfDistancesVec
	      << std::endl;
  }
  
  // Test independence of 'distance' w.r.t. order of cdfs
  horizontalDistances(valQoiRoutine_Data.qoiRv->cdf(),
		      calQoiRoutine_Data.qoiRv->cdf(),
		      *epsilonVec,
		      cdfDistancesVec);
  if (env.rank() == 0) {
    std::cout << "For epsilonVec = "    << *epsilonVec
	      << ", cdfDistancesVec (swithced order of cdfs) = " << cdfDistancesVec
	      << std::endl;
  }
  delete epsilonVec;

  // done with everything
  if (env.rank() == 0) {
    std::cout << "Finishing run of 'uqBurgers_No_Model_Uncertainty' example"
              << std::endl;
  }

  return;
}
#endif // __UQ_BURGERS_GP_MODEL_UNCERTAINTY_H__
