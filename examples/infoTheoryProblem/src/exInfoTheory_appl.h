/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef EX_INFO_THEORY_APPL_H
#define EX_INFO_THEORY_APPL_H

#include <queso/Defines.h>
#ifdef QUESO_HAS_ANN

#include <queso/StatisticalForwardProblem.h>
#include <queso/InfoTheory.h>
#include <gsl/gsl_linalg.h>
//#include <exInfoTheory_qoi.h>

static const double Pi = 3.1415926535897932;

#define length(a) ( sizeof ( a ) / sizeof ( *a ) )

//********************************************************
// The driving routine: called by main()
//********************************************************
template<class P_V,class P_M, class Q_V, class Q_M>
void 
uqAppl(const uqBaseEnvironment& env)
{
  if (env.fullRank() == 0) {
    std::cout << "Beginning run of 'exInfoTheory_example'\n"
              << std::endl;
  }


  /*********************************************
   * Entropy example
   *********************************************/
  if (env.fullRank() == 0) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "EX01: Estimate the entropy of a xD Gaussian RV" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Dim \t Analytical \t Estimated \t Rel.Abs.Err[%]" << std::endl;
  }


  for( unsigned int dim = 1; dim <= 10; ++dim )
    {

      //******************************************************
      // Ex 1: estimate the entropy of a xD Gaussian RV 
      // x = 1,2, ... ,10
      //******************************************************



      uqVectorSpace<P_V,P_M> paramSpace(env,"param_",dim,NULL);
      
      // mean vector
      P_V meanVector( paramSpace.zeroVector() );
      
      // (co)variance
      P_V varVector( paramSpace.zeroVector() );
      double prodDiag = 1.0;
      for( unsigned int i = 0; i < dim; i++ ) {
	varVector[ i ] = 3.0 + (double)i;
	prodDiag *= varVector[ i ];
      }

      // create a Gaussian RV
      uqGaussianVectorRV<P_V,P_M> gaussRV( "gaussRV_", paramSpace, meanVector, varVector );

      // plot results
      double exact_entropy = log( sqrt( pow( 2.0*Pi, (double)dim ) * prodDiag ) ) + (double)dim / 2.0;

      double est_entropy = gaussRV.estimateENT_ANN();

      if (env.fullRank() == 0) {
	std::cout << dim << " \t " << exact_entropy << " \t " << est_entropy << 
	  " \t " << fabs(exact_entropy - est_entropy) / fabs(exact_entropy) * 100.0  << std::endl;
      }

    }
  
  std::cout << "-----------------------------------------------" << std::endl;

  /*********************************************
   * Kullback-Leibler example
   *********************************************/
  if (env.fullRank() == 0) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "EX02: Estimate the Kullback-Leibler divergence" << std::endl;
    std::cout << "between two Gaussian RV" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "No.Smp. \t Analytical \t Estimated \t Rel.Abs.Err[%]" << std::endl;
  }

  // parameters
  unsigned int allN[] = {1000, 10000, 50000, 100000, 150000, 200000,
			 250000, 300000, 500000, 900000};
  unsigned int k = 1;
  double eps = UQ_INFTH_ANN_EPS;

  uqVectorSpace<P_V,P_M> paramSpace(env, "param_", 1, NULL);

  // means
  double muX = 2.0;
  P_V meanVectorX( paramSpace.zeroVector() );
  meanVectorX[ 0 ] = muX;
  double muY = 3.0;
  P_V meanVectorY( paramSpace.zeroVector() );
  meanVectorY[ 0 ] = muY;

  // covariance matrix
  double varX = pow( 1.0, 2.0 );
  P_M* covMatrixX = paramSpace.newMatrix();
  (*covMatrixX)(0,0) = varX;
  double varY = pow( 1.0, 2.0 );
  P_M* covMatrixY = paramSpace.newMatrix();
  (*covMatrixY)(0,0) = varY;

  // create Gaussian RV
  uqGaussianVectorRV<P_V,P_M> xRV( "xRV_", paramSpace, meanVectorX, *covMatrixX );
  uqGaussianVectorRV<P_V,P_M> yRV( "yRV_", paramSpace, meanVectorY, *covMatrixY );

  double exact_kl = 0.5 * ( log(varY/varX) + varX/varY + pow(muY-muX,2.0)/varY - 1.0 );
  
  unsigned int xDimSel[1] = { 0 };
  unsigned int yDimSel[1] = { 0 }; 

  for( unsigned int indN = 0; indN < length(allN); ++indN )
    {
      unsigned int xN = allN[ indN ]; 
      unsigned int yN = allN[ indN ];
      double est_kl = estimateKL_ANN( xRV, yRV, xDimSel, 1, yDimSel, 1, xN, yN, k, eps );

      if (env.fullRank() == 0) {
	std::cout << xN << " \t\t " << exact_kl << " \t " << est_kl << 
	  " \t " << fabs(exact_kl - est_kl) / fabs(exact_kl) * 100.0  << std::endl;
      }
    }

  std::cout << "-----------------------------------------------" << std::endl;

  /*********************************************
   * Cross-entropy example
   *********************************************/
  if (env.fullRank() == 0) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "EX03: Estimate the Cross-Entropy" << std::endl;
    std::cout << "between two Gaussian RV" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "No.Smp. \t Analytical \t Estimated \t Rel.Abs.Err[%]" << std::endl;

  }

  double exact_ce = 0.5 * ( log(2.0*Pi) + log(varY) + pow(muY,2.0)/varY + varX/varY +
			    pow(muX,2.0)/varY - 2.0*muX*muY/varY);
  
  for( unsigned int indN = 0; indN < length(allN); ++indN )
    {
      unsigned int xN = allN[ indN ]; 
      unsigned int yN = allN[ indN ];
      double est_ce = estimateCE_ANN( xRV, yRV, xDimSel, 1, yDimSel, 1, xN, yN, k, eps );

      if (env.fullRank() == 0) {
	std::cout << xN << " \t\t " << exact_ce << " \t " << est_ce << 
	  " \t " << fabs(exact_ce - est_ce) / fabs(exact_ce) * 100.0  << std::endl;
      }
    }

  std::cout << "-----------------------------------------------" << std::endl;


  /*********************************************
   * Mutual information example
   *********************************************/
  if (env.fullRank() == 0) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "EX04: Estimate the Mutual information" << std::endl;
    std::cout << "between two Gaussian RV" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "No.Smp. \t Analytical \t Estimated \t Rel.Abs.Err[%]" << std::endl;

  }

  unsigned int allN_MI[] = {1000, 10000, 50000, 100000};
  k = UQ_INFTH_ANN_KNN;
  eps = UQ_INFTH_ANN_EPS;

  unsigned int dim = 2;
  uqVectorSpace<P_V,P_M> paramSpace_MI(env, "param_", dim, NULL);

  // means
  P_V meanVector( paramSpace_MI.zeroVector() );
  meanVector[ 0 ] = 2.0;
  meanVector[ 1 ] = 4.0;

  // covariance matrix
  double rho = 0.8;
  double std1 = 1.0;
  double std2 = 1.0;
  P_M* covMatrix = paramSpace_MI.newMatrix();
  (*covMatrix)(0,0) = std1*std1; (*covMatrix)(0,1) = rho*std1*std2;
  (*covMatrix)(1,0) = rho*std1*std2; (*covMatrix)(1,1) = std2*std2;

  // create a Gaussian RV
  uqGaussianVectorRV<P_V,P_M> gaussRV( "gaussRV_", paramSpace_MI, meanVector, *covMatrix );

  // get the determinant of the matrix
  int signum;
  gsl_permutation* p = gsl_permutation_alloc( dim );
  gsl_matrix* mat = covMatrix->data();
  gsl_linalg_LU_decomp ( mat, p, &signum );
  double matDet = gsl_linalg_LU_det( mat, signum );

  double exact_mutinfo = - 0.5 * log( matDet );

  unsigned int xDimSel_MI[1] = { 0 };
  unsigned int yDimSel_MI[1] = { 1 }; 

  for( unsigned int indN = 0; indN < length(allN_MI); ++indN )
    {
      unsigned int N = allN_MI[ indN ]; 
      double est_mutinfo = estimateMI_ANN( gaussRV, xDimSel_MI, 1, yDimSel_MI, 1, k, N, eps );

      if (env.fullRank() == 0) {
	std::cout << N << " \t\t " << exact_mutinfo << " \t " << est_mutinfo <<  " \t " <<
	  fabs(exact_mutinfo - est_mutinfo) / fabs(exact_mutinfo) * 100.0 << std::endl;
      }
    }


  return;
}

#endif // QUESO_HAS_ANN

#endif // EX_INFO_THEORY_APPL_H
