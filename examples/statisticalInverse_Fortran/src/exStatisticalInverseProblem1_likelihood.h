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
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __EX_STATISTICAL_INVERSE_PROBLEM_1_LIKELIHOOD_H__
#define __EX_STATISTICAL_INVERSE_PROBLEM_1_LIKELIHOOD_H__

#include<cstdio>
#include<stdlib.h>
#include<iostream>

#define FORTRAN_WRAPPER_VERSION

// Prototype C function for passing likelihood function to Fortran

extern "C" double f_likelihood(int num_params, const double *parameter_values, 
			       const double *parameter_means, const double *matrix);

//-------------------------------------------------------------
// The likelihood routine: provided by user and called by QUESO
//-------------------------------------------------------------

template<class P_V,class P_M>
struct
likelihoodRoutine_DataType
{
  const P_V* paramMeans;
  const P_M* matrix;
  bool       applyMatrixInvert;
};

template<class P_V,class P_M>
double
likelihoodRoutine(
  const P_V&  paramValues,
  const P_V*  paramDirection,
  const void* functionDataPtr,
  P_V*        gradVector,
  P_M*        hessianMatrix,
  P_V*        hessianEffect)
{

  const P_V& paramMeans        = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->paramMeans;
  const P_M& matrix            = *((likelihoodRoutine_DataType<P_V,P_M> *) functionDataPtr)->matrix;

#ifdef FORTRAN_WRAPPER_VERSION

  // grab data in a form that we can pass to fortran

  const double *params      = gsl_vector_const_ptr(paramValues.data(),0);
  const double *param_means = gsl_vector_const_ptr(paramMeans.data(),0);

  static double *matrix_to_fortran;
  static bool first_entry=true;

  // Initialization to perform on first entry only

  if(first_entry)
    {
      matrix_to_fortran = (double *)calloc(matrix.numRowsGlobal()*matrix.numCols(),sizeof(double));

      for(int row=0;row<matrix.numRowsGlobal();row++)
	for(int col=0;col<matrix.numCols();col++)
	  {
	    matrix_to_fortran[col+row*matrix.numCols()] = matrix(row,col);
	  }

      first_entry=false;
    }

  // Get likelihood from Fortran function

  return(f_likelihood(paramValues.sizeGlobal(),params,param_means,matrix_to_fortran));

#else

  double result = 0.;

  P_V diffVec(paramValues - paramMeans);
  result = scalarProduct(diffVec, matrix * diffVec);

  return(result);
#endif

}
#endif // __EX_STATISTICAL_INVERSE_PROBLEM_1_LIKELIHOOD_H__
