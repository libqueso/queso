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
  uqEnvironmentClass* env = new uqEnvironmentClass(argc,argv);

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
calibrationCovarianceChol(const int N, const double *xLoc, const double sig2, const double ellx, 
			  double *cholK)
{
  int ierr;
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
