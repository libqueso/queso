/* uq/examples/queso/burgers/no_model_uncertainty/uqBurgers_GP_Model_Uncertainty_gsl.C
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

#include <uqBurgers_GP_Model_Uncertainty.h>
#include <uqGslMatrix.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environment
  //************************************************
  MPI_Init(&argc,&argv);
  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(argc,argv);

  //************************************************
  // Call application
  //************************************************
  uqAppl<uqGslVectorClass, // type for parameter vectors
         uqGslMatrixClass, // type for parameter matrices
         uqGslVectorClass, // type for qoi vectors
         uqGslMatrixClass  // type for qoi matrices
        >(*env);

  //************************************************
  // Finalize environment
  //************************************************
  delete env;
  MPI_Finalize();

  return 0;
}


// Computes covariance matrix corresponding to the 
// squared exponential covariance function, and
// uses to compute log(likelihood)
void
negTwoLogLikelihood(const int N, const double *cholK, const double *misfit, double *result)
{
  int ierr;
  double mt_iK_m, hlndet;
  double alfa[N];

  // compute K\misfit using previously computed chol decomposition
  gsl_matrix_const_view gslK = gsl_matrix_const_view_array(cholK, N, N);
  gsl_vector_const_view gslMis = gsl_vector_const_view_array(misfit, N);
  gsl_vector_view gslAlfa = gsl_vector_view_array(alfa, N);

  ierr = gsl_linalg_cholesky_solve(&gslK.matrix, &gslMis.vector, &gslAlfa.vector);
  if( ierr == GSL_EDOM ){
    printf("ERROR: MATRIX IS NOT POSITIVE DEFINITE!\n");
    return;
  }

  // misfit'*(K\misfit)
  mt_iK_m = 0.0;
  for( int ii=0; ii<N; ii++ ) mt_iK_m += misfit[ii]*alfa[ii];
  
  // Compute the 0.5*log(det(K))
  hlndet = log(cholK[0]);
  for( int ii=1; ii<N; ii++ ) hlndet += log(cholK[N*ii+ii]);

  // Compute log(likelihood)
  *result = mt_iK_m + 2.0*hlndet + N*log(2.0*PI);

  return;
}


// Computes Cholesky decomposition of GP prior covariance
// required for calibration phase 
int
calibrationCovarianceChol(const double sig2, const double ellx, inverseProblem_DataClass *cal_data)
{
  int ierr, N=cal_data->nDataPoints;
  const double *xLoc = cal_data->dataLocations;
  double *cholK = cal_data->cholK;
  const double noise = 1e-10;

  // Compute covariance matrix
  for( int ii=0; ii<N; ii++ ){
    for( int jj=0; jj<N; jj++ ){
      double dxol = (xLoc[ii] - xLoc[jj])/ellx;
      cholK[N*ii+jj] = sig2*exp(-0.5*dxol*dxol); //K(ii, jj) = sig2*exp(-0.5*dxol*dxol);
    }
    cholK[N*ii+ii] += noise; //K(ii, ii) += noise; // add noise to diagonal to eliminate ill-conditioning
  }

  // Compute Cholesky decomposition
  // BELOW HERE, K HAS BEEN OVERWRITTEN WITH chol(K)
  gsl_matrix_view gslCholK = gsl_matrix_view_array(cholK, N, N);
  ierr = gsl_linalg_cholesky_decomp(&gslCholK.matrix);
  if( ierr == GSL_EDOM ){
    printf("ERROR: MATRIX IS NOT POSITIVE DEFINITE!\n");
    return ierr;
  }

  return 0;
}

// Computes GP prior covariance between calibration and validation data points
int
computePriorCovarianceOffDiag(const int N1, const double *xLoc1, const double Re1,
			      const int N2, const double *xLoc2, const double Re2,
			      const double sig2, const double ellx, const double ellRe,
			      double *K21)
{
  const double logRe1 = log10(Re1);
  const double logRe2 = log10(Re2);
  const double dlogReol = (logRe2 - logRe1)/ellRe;
  const double dlogReol2 = dlogReol*dlogReol;
  double dxol;

  // Compute covariance matrix
  for( int ii=0; ii<N2; ii++ ){
    for( int jj=0; jj<N1; jj++ ){
      dxol = (xLoc2[ii] - xLoc1[jj])/ellx;
      K21[N1*ii+jj] = sig2*exp( -0.5*(dxol*dxol + dlogReol2) );
    }
  }

  return 0;
}


// Compute GP mean for validation phase (i.e. posterior mean from calibration)
int
validationMean(const double sig2, const double ellx, const double ellRe, const double kappa,
	       likelihoodRoutine_DataClass *data, double *mean)
{

  int ierr;
  int N1 = data->calibrationData->nDataPoints, N2 = data->validationData->nDataPoints;
  double cal_umodel[N1];
  const double *cal_udata = data->calibrationData->dataValues;
  
  double Re1 = data->calibrationData->Re, Re2 = data->validationData->Re;
  double *xLoc1 = data->calibrationData->dataLocations, *xLoc2 = data->validationData->dataLocations;
  
//   quadBasis *pQB = data->pQB;
//   gsl_vector *U = data->U;
  
//   // evaluate misfit for stage 1 scenario
//   ierr = solveForStateAtXLocations(xLoc1, N1, 1.0/Re1, kappa, U, pQB, cal_umodel);
//   if( ierr != 0 ){
//     printf("WARNING: Burgers solver was not successful!\n");
//     printf("         Results are probably meaningless.\n");
//     fflush(stdout);
//   }

  burgersQuesoInterface& calInterface = data->calibrationBurgersInterface;
  ierr = calInterface.solveForStateAtDataLocations(1.0/Re1, kappa, cal_umodel);
  UQ_FATAL_TEST_MACRO( (ierr!=0), UQ_UNAVAILABLE_RANK,
		       "uqAppl(), in uqBurgers_No_Model_Uncertainty (likelihoodRoutine)",
		       "failed solving Burgers' eqn" );
    
  double cal_misfit[N1];
  for( int ii=0; ii<N1; ii++ ) cal_misfit[ii] = (cal_udata[ii] - cal_umodel[ii]);


  double *cal_cholK = data->calibrationData->cholK;

  double K21[N2*N1];
  double temp[N1];

  // Compute stage 1 off diagonal
  ierr = computePriorCovarianceOffDiag(N1, xLoc1, Re1, N2, xLoc2, Re2, sig2, ellx, ellRe, K21);
  if( ierr != 0 ) return ierr;

  // temp = K\s1_misfit
  gsl_matrix_const_view gslK = gsl_matrix_const_view_array(cal_cholK, N1, N1);
  gsl_vector_const_view gslMis = gsl_vector_const_view_array(cal_misfit, N1);
  gsl_vector_view gslTemp = gsl_vector_view_array(temp, N1);

  ierr = gsl_linalg_cholesky_solve(&gslK.matrix, &gslMis.vector, &gslTemp.vector);
  if( ierr == GSL_EDOM ){
    printf("ERROR: MATRIX IS NOT POSITIVE DEFINITE!\n");
    return ierr;
  }

  // s1_mean = K21*(K\s1_misfit)
  gsl_matrix_const_view gslK21 = gsl_matrix_const_view_array(K21, N2, N1);
  gsl_vector_view gslMean = gsl_vector_view_array(mean, N2);

  ierr = gsl_blas_dgemv(CblasNoTrans, 1.0, &gslK21.matrix, &gslTemp.vector, 0.0, &gslMean.vector);
  if( ierr != 0 ) return ierr;

  return 0;
}


// Computes Cholesky decomposition of GP covariance
// required for stage 2 calibration
int
validationCovarianceChol(const double sig2, const double ellx, const double ellRe,
			 const inverseProblem_DataClass *cal_data, inverseProblem_DataClass *val_data)
{
  int ierr;
  int N1 = cal_data->nDataPoints, N2 = val_data->nDataPoints;
  
  double Re1 = cal_data->Re, Re2 = val_data->Re;
  double *xLoc1 = cal_data->dataLocations, *xLoc2 = val_data->dataLocations;

  double *cal_cholK = cal_data->cholK, *val_cholK = val_data->cholK;

  double K21[N2*N1];
  double temp[N1*N2];

  // Compute stage 1 off diagonal
  ierr = computePriorCovarianceOffDiag(N1, xLoc1, Re1, N2, xLoc2, Re2, sig2, ellx, ellRe, K21);
  if( ierr != 0 ) return ierr;

  // temp = transpose(K21)
  for( int ii=0; ii<N1; ii++ ){
    for( int jj=0; jj<N2; jj++ ){
      temp[N2*ii+jj] = K21[N1*jj+ii];
    }
  }

  // temp = (K11 + M)^{-1}*K12
  gsl_matrix_const_view gslCalCholK = gsl_matrix_const_view_array(cal_cholK, N1, N1);
  gsl_matrix_view gslTemp = gsl_matrix_view_array(temp, N1, N2);

  ierr = gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, &gslCalCholK.matrix, &gslTemp.matrix);
  if( ierr != 0 ) return ierr;

  ierr = gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &gslCalCholK.matrix, &gslTemp.matrix);
  if( ierr != 0 ) return ierr;
  
  // Compute stage 2 covariace according to prior
  for( int ii=0; ii<N2; ii++ ){
    for( int jj=0; jj<N2; jj++ ){
      double dxol = (xLoc2[ii] - xLoc2[jj])/ellx;
      val_cholK[N2*ii+jj] = sig2*exp(-0.5*dxol*dxol); //K(ii, jj) = sig2*exp(-0.5*dxol*dxol);
    }
  }

  // C = K22 - K21*(K11 + noise_variance*I)^{-1}*K12
  gsl_matrix_view gslValCholK = gsl_matrix_view_array(val_cholK, N2, N2);
  gsl_matrix_view gslK21 = gsl_matrix_view_array(K21, N2, N1);

  ierr = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &gslK21.matrix, &gslTemp.matrix, 1.0, &gslValCholK.matrix);
  if( ierr != 0 ) return ierr;

  // Add noise to diagonal to ensure Cholesky decomposition is successful
  for( int ii=0; ii<N2; ii++ ){
    val_cholK[N2*ii+ii] += 1e-10;
  }
      

  // Cholesky decomposition
  ierr = gsl_linalg_cholesky_decomp(&gslValCholK.matrix);
  if( ierr == GSL_EDOM ){
    printf("ERROR: MATRIX IS NOT POSITIVE DEFINITE!\n");
    return ierr;
  }

  return 0;
}


int
computeStage1PropagationMean(const int N1, const double *xLoc1, const double Re1,
			     const int N3, const double *xLoc3, const double Re3,
			     const double sig2, const double ellx, const double ellRe,
			     const double *s1_misfit, const double *cholK, double *mean)
{

  // variance = K33_xx - K13_x^T*(M+K11)^-1*K13_x
  double K13[N1], K13_x[N1];
  double dxol, dlogRe;

  dlogRe = (log10(Re3) - log10(Re1))/ellRe;

  for( int ii=0; ii<N1; ii++ ){
    dxol = (xLoc1[ii] - xLoc3[0])/ellx;
    K13[ii] = sig2*exp(-0.5*(dxol*dxol + dlogRe*dlogRe) );
    K13_x[ii] = (dxol/ellx)*K13[ii];
  }
  
  gsl_matrix_const_view gslCholK = gsl_matrix_const_view_array(cholK, N1, N1);
  gsl_vector_const_view gslMis = gsl_vector_const_view_array(s1_misfit, N1);

  double tmp[N1];
  gsl_vector_view gslTmp = gsl_vector_view_array(tmp, N1);

  int ierr = gsl_linalg_cholesky_solve(&gslCholK.matrix, &gslMis.vector, &gslTmp.vector);
  if( ierr != 0 ) return ierr;

  *mean=0.0;
  for( int ii=0; ii<N1; ii++ ) *mean += K13_x[ii]*tmp[ii];

  return 0;
}

int
computeStage1PropagationVariance(const int N1, const double *xLoc1, const double Re1,
				 const int N3, const double *xLoc3, const double Re3,
				 const double sig2, const double ellx, const double ellRe,
				 const double *cholK, double *variance)
{

  // variance = K33_xx - K13_x^T*(M+K11)^-1*K13_x
  double K33 = sig2;
  double K33_xx = K33/(ellx*ellx);

  double K13[N1], K13_x[N1];
  double dxol, dlogRe;

  dlogRe = (log10(Re3) - log10(Re1))/ellRe;

  for( int ii=0; ii<N1; ii++ ){
    dxol = (xLoc1[ii] - xLoc3[0])/ellx;
    K13[ii] = sig2*exp(-0.5*(dxol*dxol + dlogRe*dlogRe) );
    K13_x[ii] = (dxol/ellx)*K13[ii];
  }
  
  gsl_matrix_const_view gslCholK = gsl_matrix_const_view_array(cholK, N1, N1);
  gsl_vector_const_view gslK13_x = gsl_vector_const_view_array(K13_x, N1);

  double tmp[N1];
  gsl_vector_view gslTmp = gsl_vector_view_array(tmp, N1);

  int ierr = gsl_linalg_cholesky_solve(&gslCholK.matrix, &gslK13_x.vector, &gslTmp.vector);
  if( ierr != 0 ) return ierr;

  double sum=0.0;
  for( int ii=0; ii<N1; ii++ ) sum += K13_x[ii]*tmp[ii];

  (*variance) = K33_xx - sum;

  return 0;
}



int
computeStage2PropagationMean(const int N1, const double *xLoc1, const double Re1,
			     const int N2, const double *xLoc2, const double Re2,
			     const int N3, const double *xLoc3, const double Re3,
			     const double sig2, const double ellx, const double ellRe,
			     const double *cholK1, const double *cholK2, double *s2_misfit, double mean1,
			     double *mean2)
{
  int ierr;

  // compute K12
  double K21[N2*N1];
  ierr = computePriorCovarianceOffDiag(N1, xLoc1, Re1, N2, xLoc2, Re2, sig2, ellx, ellRe, K21);
  if( ierr != 0 ) return ierr;

  // temp = transpose(K21)
  double temp[N1*N2];
  for( int ii=0; ii<N1; ii++ ){
    for( int jj=0; jj<N2; jj++ ){
      temp[N2*ii+jj] = K21[N1*jj+ii];
    }
  }

  // temp = (K11 + M)^{-1}*K12
  gsl_matrix_const_view gslCholK1 = gsl_matrix_const_view_array(cholK1, N1, N1);
  gsl_matrix_view gslTemp = gsl_matrix_view_array(temp, N1, N2);

  ierr = gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, &gslCholK1.matrix, &gslTemp.matrix);
  if( ierr != 0 ) return ierr;

  ierr = gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &gslCholK1.matrix, &gslTemp.matrix);
  if( ierr != 0 ) return ierr;

  // compute K31_x
  double K13[N1], K13_x[N1];
  double dxol, dlogRe;

  dlogRe = (log10(Re3) - log10(Re1))/ellRe;

  for( int ii=0; ii<N1; ii++ ){
    dxol = (xLoc1[ii] - xLoc3[0])/ellx;
    K13[ii] = sig2*exp(-0.5*(dxol*dxol + dlogRe*dlogRe) );
    K13_x[ii] = (dxol/ellx)*K13[ii];
  }

  // temp2 = K31_x*(K11 + M)^{-1}*K12
  double temp2[N2];
  for( int ii=0; ii<N2; ii++ ){
    temp2[ii] = 0.0;
    for( int jj=0; jj<N1; jj++ ){
      temp2[ii] += K13_x[jj]*temp[N2*jj+ii];
    }
  }

  // Compute K32_x
  double K23[N2], K23_x[N2];

  dlogRe = (log10(Re3) - log10(Re2))/ellRe;

  for( int ii=0; ii<N2; ii++ ){
    dxol = (xLoc2[ii] - xLoc3[0])/ellx;
    K23[ii] = sig2*exp(-0.5*(dxol*dxol + dlogRe*dlogRe) );
    K23_x[ii] = (dxol/ellx)*K23[ii];
  }

  // Compute C32
  double C32[N2];
  for( int ii=0; ii<N2; ii++ ) C32[ii] = K23_x[ii] - temp2[ii];


  gsl_matrix_const_view gslCholK2 = gsl_matrix_const_view_array(cholK2, N2, N2);
  gsl_vector_const_view gslMis = gsl_vector_const_view_array(s2_misfit, N2);

  double tmp[N2];
  gsl_vector_view gslTmp = gsl_vector_view_array(tmp, N2);

  ierr = gsl_linalg_cholesky_solve(&gslCholK2.matrix, &gslMis.vector, &gslTmp.vector);
  if( ierr != 0 ) return ierr;

  *mean2=0.0;
  for( int ii=0; ii<N2; ii++ ) *mean2 += C32[ii]*tmp[ii];
  *mean2 += mean1;

  return 0;
}


int
computeStage2PropagationVariance(const int N1, const double *xLoc1, const double Re1,
				 const int N2, const double *xLoc2, const double Re2,
				 const int N3, const double *xLoc3, const double Re3,
				 const double sig2, const double ellx, const double ellRe,
				 const double *cholK1, const double *cholK2, const double variance1,
				 double *variance2)
{

  int ierr;

  // compute K12
  double K21[N2*N1];
  ierr = computePriorCovarianceOffDiag(N1, xLoc1, Re1, N2, xLoc2, Re2, sig2, ellx, ellRe, K21);
  if( ierr != 0 ) return ierr;

  // temp = transpose(K21) = K12
  double temp[N1*N2];
  for( int ii=0; ii<N1; ii++ ){
    for( int jj=0; jj<N2; jj++ ){
      temp[N2*ii+jj] = K21[N1*jj+ii];
    }
  }

  // temp = (K11 + M)^{-1}*K12
  gsl_matrix_const_view gslCholK1 = gsl_matrix_const_view_array(cholK1, N1, N1);
  gsl_matrix_view gslTemp = gsl_matrix_view_array(temp, N1, N2);

  ierr = gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, &gslCholK1.matrix, &gslTemp.matrix);
  if( ierr != 0 ) return ierr;

  ierr = gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &gslCholK1.matrix, &gslTemp.matrix);
  if( ierr != 0 ) return ierr;

  // compute K31_x
  double K13[N1], K13_x[N1];
  double dxol, dlogRe;

  dlogRe = (log10(Re3) - log10(Re1))/ellRe;

  for( int ii=0; ii<N1; ii++ ){
    dxol = (xLoc1[ii] - xLoc3[0])/ellx;
    K13[ii] = sig2*exp(-0.5*(dxol*dxol + dlogRe*dlogRe) );
    K13_x[ii] = (dxol/ellx)*K13[ii];
    //printf("xx1[%d] = %.6E, xx3 = %.6E, Re1 = %.6E, Re3 = %.6E\n", ii, xLoc1[ii], xLoc3[0], Re1, Re3);
    //printf("K13[%d] = %.6E, K13_x[%d] = %.6E\n", ii, K13[ii], ii, K13_x[ii]); fflush(stdout);
    //K13_x[ii] = -(dxol/ellx)*K13[ii];
  }

  // temp2 = K31_x*(K11 + M)^{-1}*K12
  double temp2[N2];
  for( int ii=0; ii<N2; ii++ ){
    temp2[ii] = 0.0;
    for( int jj=0; jj<N1; jj++ ){
      temp2[ii] += K13_x[jj]*temp[N2*jj+ii];
    }
  }

  // Compute K32_x
  double K23[N2], K23_x[N2];

  dlogRe = (log10(Re3) - log10(Re2))/ellRe;

  for( int ii=0; ii<N2; ii++ ){
    dxol = (xLoc2[ii] - xLoc3[0])/ellx;
    K23[ii] = sig2*exp(-0.5*(dxol*dxol + dlogRe*dlogRe) );
    K23_x[ii] = (dxol/ellx)*K23[ii];
    //K23_x[ii] = -(dxol/ellx)*K23[ii];
  }

  // Compute C32
  double C32[N2];
  for( int ii=0; ii<N2; ii++ ){
    printf("K23_x[%d] = %.6E, temp2[%d] = %.6E\n", ii, K23_x[ii], ii, temp2[ii]); fflush(stdout);
    C32[ii] = K23_x[ii] - temp2[ii];
  }


  gsl_matrix_const_view gslCholK2 = gsl_matrix_const_view_array(cholK2, N2, N2);
  gsl_vector_const_view gslC32 = gsl_vector_const_view_array(C32, N2);

  double tmp[N2];
  gsl_vector_view gslTmp = gsl_vector_view_array(tmp, N2);

  ierr = gsl_linalg_cholesky_solve(&gslCholK2.matrix, &gslC32.vector, &gslTmp.vector);
  if( ierr != 0 ) return ierr;

  double sum=0.0;
  for( int ii=0; ii<N2; ii++ ){
    printf("C32[%d] = %.6E, tmp[%d] = %.6E\n", ii, C32[ii], ii, tmp[ii]); fflush(stdout);
    sum += C32[ii]*tmp[ii];
  }


  //printf("variance1 = %.6E, sum = %.6E\n", variance1, sum); fflush(stdout);

  (*variance2) = variance1 - sum;

  return 0;
}
