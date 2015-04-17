//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef EX_STATISTICAL_INVERSE_PROBLEM_1_LIKELIHOOD_H
#define EX_STATISTICAL_INVERSE_PROBLEM_1_LIKELIHOOD_H

#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <queso/Environment.h>

#define FORTRAN_WRAPPER_VERSION

// Prototype C function for passing likelihood function to Fortran

extern "C" double f_likelihood(int num_params, const double *parameter_values,
			       const double *parameter_means, const double *matrix, int subid,
			       MPI_Fint subcomm);

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

  static int num_procs, num_local;
  static int subId;
  static MPI_Comm Local_MPI_Comm;

  if(first_entry)
    {
      matrix_to_fortran = (double *)calloc(matrix.numRowsGlobal()*matrix.numCols(),sizeof(double));

      for(unsigned int row=0;row<matrix.numRowsGlobal();row++)
	for(unsigned int col=0;col<matrix.numCols();col++)
	  {
	    matrix_to_fortran[col+row*matrix.numCols()] = matrix(row,col);
	  }

      // determine local MPI environment info

      const QUESO::BaseEnvironment &env = paramMeans.env();
      subId                             = (int)env.subId();
      Local_MPI_Comm                    = env.subComm().Comm();
      int num_global_procs;

      MPI_Comm_size (Local_MPI_Comm, &num_procs);
      MPI_Comm_rank (Local_MPI_Comm, &num_local);

      if(paramMeans.env().worldRank() == 0)
	{
	  MPI_Comm_size (MPI_COMM_WORLD, &num_global_procs);
	  printf("\nLikelihood parallel environment (C++): \n");
	  printf("--> Total # of MPI processors available    = %i\n",  num_global_procs);
	  printf("--> Total # of SubEnvironments requested   = %i\n",  paramMeans.env().numSubEnvironments());
	  printf("--> Number of MPI tasks per SubEnvironment = %i\n\n",num_procs);
	}

      first_entry=false;
    }

  // Transfer program to Fortran likelihood function - note that rank 0
  // for each subcommunicator is responsible for computing the final
  // likelihood value.

  double result = f_likelihood(paramValues.sizeGlobal(),params,param_means,matrix_to_fortran,subId,
			       MPI_Comm_c2f(Local_MPI_Comm));

  return(result);

#else

  double result = 0.;

  P_V diffVec(paramValues - paramMeans);
  result = scalarProduct(diffVec, matrix * diffVec);

  return(result);
#endif

}
#endif // EX_STATISTICAL_INVERSE_PROBLEM_1_LIKELIHOOD_H
