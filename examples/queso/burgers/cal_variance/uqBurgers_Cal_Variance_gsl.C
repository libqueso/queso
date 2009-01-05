/* uq/examples/queso/burgers/embedded_model_uncertainty/uqBurgers_Cal_Variance_gsl.C
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

#include <iostream>

#include <uqBurgers_Cal_Variance.h>
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


int
computeBaseModelCovariance(const int nScenario, const double *Re, const int *nXloc, double **xx,
			   double **K)
{
  // compute total number of data points
  int nTot = 0;
  for( int ii=0; ii<nScenario; ii++ ) nTot += nXloc[ii];

  // allocate covariance
  *K = new double[nTot*nTot];

  // compute covariance
  int row = 0;
  int col = 0;

  double iRe, jRe, dReol;
  double ix, jx, dxol;

  for( int iScen=0; iScen<nScenario; iScen++ ){
    iRe = Re[iScen];
    for( int iPt=0; iPt<nXloc[iScen]; iPt++ ){
      ix = xx[iScen][iPt];

      col = 0;
      for( int jScen=0; jScen<nScenario; jScen++ ){
	jRe = Re[jScen];
	
	dReol = (log10(iRe) - log10(jRe))/0.5;

	for( int jPt=0; jPt<nXloc[jScen]; jPt++ ){
	  jx = xx[jScen][jPt];

	  dxol = (ix - jx)/0.1;

	  //std::cout << "(ix, jx) = (" << ix << ", " << jx << ")\n" << std::flush;


	  (*K)[nTot*row + col] = exp(-0.5*(dxol*dxol + dReol*dReol));
// 	  std::cout << "K(" << row << ", " << col << ") = " << 
// 	    std::setprecision(15) << std::scientific << (*K)[nTot*row + col] << "\n" << std::flush;


	  col++;
	  //std::cout << "col = " << col << "\n" << std::flush;
	}
      }
      
      row++;
      //std::cout << "row = " << row << "\n" << std::flush;
    }
  }

  //std::cout << "(row, col) = (" << row << ", " << col << ")\n" << std::flush;

  return 0;
}

int
computeLikelihood(likelihoodRoutine_DataClass<uqGslVectorClass, uqGslMatrixClass>* functionDataPtr, 
		  double modelVar, double expVar, double* misfit, double& nTwoLogLikelihood)
{
  // get data
  const int nScenario = ((likelihoodRoutine_DataClass<uqGslVectorClass, uqGslMatrixClass> *) functionDataPtr)->nScenario;
  const int *nXloc = ((likelihoodRoutine_DataClass<uqGslVectorClass, uqGslMatrixClass> *) functionDataPtr)->nXloc;
  double **xx = ((likelihoodRoutine_DataClass<uqGslVectorClass, uqGslMatrixClass> *) functionDataPtr)->xx;
  const double *Re = ((likelihoodRoutine_DataClass<uqGslVectorClass, uqGslMatrixClass> *) functionDataPtr)->Re;

  const double myPI = 3.14159265358979323846;

  int ierr;
  double *K;

  // compute total number of data points
  int nTot = 0;
  for( int ii=0; ii<nScenario; ii++ ) nTot += nXloc[ii];

  // Compute basic model covariance matrix
  ierr = computeBaseModelCovariance(nScenario, Re, nXloc, xx, &K);
  if( ierr != 0 ) return ierr;

  // Multiply basic model covariance by model variance
  for( int ii=0; ii<nTot*nTot; ii++ ) K[ii] *= modelVar;

  // Add experimental variance to the diagonal
  for( int ii=0; ii<nTot; ii++ ) K[nTot*ii + ii] += expVar;


  // Compute Cholesky decomposition of covariance
  gsl_matrix_view gslK = gsl_matrix_view_array(K, nTot, nTot);

  ierr = gsl_linalg_cholesky_decomp(&gslK.matrix);
  if( ierr == GSL_EDOM ){
    printf("ERROR: MATRIX IS NOT POSITIVE DEFINITE!\n");
    return ierr;
  }

  // Compute misfit portion of likelihood
  gsl_vector *tmpVec = gsl_vector_calloc(nTot);
  gsl_vector_view gslMisfit = gsl_vector_view_array(misfit, nTot);

  // FIX ME: FINISH misfit'*(K\misfit)
  ierr = gsl_linalg_cholesky_solve(&gslK.matrix, &gslMisfit.vector, tmpVec);
  if( ierr == GSL_EDOM ){
    printf("ERROR: MATRIX IS NOT POSITIVE DEFINITE!\n");
    return ierr;
  }
  
  double mt_iK_m = 0.0;
  gsl_blas_ddot(&gslMisfit.vector, tmpVec, &mt_iK_m);

  // Compute -log(det(K))
  double logDetK = 0.0;
  for( int ii=0; ii<nTot; ii++ ) logDetK += 2.0*log(K[nTot*ii+ii]);
  
  // Compute -2*log(likelihood)... using matlab notation:
  // -2*log(likelihood) = misfit'*(K\misfit) + log(det(K)) + N*log(2*pi)
  nTwoLogLikelihood = mt_iK_m + logDetK + nTot*log(2.0*myPI);

  // clean up
  delete[] K;

  return 0;
}


