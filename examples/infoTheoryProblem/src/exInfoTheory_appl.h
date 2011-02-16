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
 *
 * $Id: exInfoTheory_appl.h 14596 2011-02-03 15:10:42Z gabriel $
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __EX_INFO_THEORY_APPL_H__
#define __EX_INFO_THEORY_APPL_H__

#include <uqDefines.h>
#ifdef QUESO_HAS_ANN

#include <uqStatisticalForwardProblem.h>
//#include <exInfoTheory_qoi.h>

static const double Pi = 3.1415926535897932;

//********************************************************
// The driving routine: called by main()
//********************************************************
template<class P_V,class P_M, class Q_V, class Q_M>
void 
uqAppl(const uqBaseEnvironmentClass& env)
{
  if (env.fullRank() == 0) {
    std::cout << "Beginning run of 'exInfoTheory_example'\n"
              << std::endl;
  }

  if (env.fullRank() == 0) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "EX01: Estimate the entropy of a xD Gaussian RV" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Dim \t Analytical \t Estimated \t Abs.Diff." << std::endl;
  }


  for( unsigned int dim = 1; dim <= 10; ++dim )
    {

      //******************************************************
      // Ex 1: estimate the entropy of a xD Gaussian RV 
      // x = 1,2, ... ,10
      //******************************************************

      uqVectorSpaceClass<P_V,P_M> paramSpace(env,"param_",dim,NULL);
      
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
      uqGaussianVectorRVClass<P_V,P_M> gaussRV( "gaussRV_", paramSpace, meanVector, varVector );

      // plot results
      double exact_entropy = log( sqrt( pow( 2.0*Pi, (double)dim ) * prodDiag ) ) + (double)dim / 2.0;

      double est_entropy = gaussRV.estimateENT_ANN();

      if (env.fullRank() == 0) {
	std::cout << dim << " \t " << exact_entropy << " \t " << est_entropy << 
	  " \t " << fabs( exact_entropy - est_entropy )  << std::endl;
      }

    }
  
  std::cout << "-----------------------------------------------" << std::endl;

  return;
}
#endif // QUESO_HAS_ANN

#endif // __EX_INFO_THEORY_APPL_H__
