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
// examples/validationCycle/src/exTgaValidationCycle_gsl.C.  It just
// generates smaller chains and so is faster and more suitable to
// serve as an automatic regression test. Indeed, it was set by Karl
// as a first example to the QUESO team on how to run regression tests
// with buildbot and how to setup them in the configure & make package
// of QUESO.
//
//--------------------------------------------------------------------------


#include <TgaValidationCycle_appl.h>
#include <queso/GslMatrix.h>

int main(int argc, char* argv[])
{
  //************************************************
  // Initialize environment
  //************************************************
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  UQ_FATAL_TEST_MACRO(argc != 2,
                      QUESO::UQ_UNAVAILABLE_RANK,
                      "main()",
                      "input file must be specified in command line as argv[1], just after executable argv[0]");
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(argv[1],"",NULL);
#endif

  //************************************************
  // Run application
  //************************************************
  uqAppl<QUESO::GslVector, // type for parameter vectors
         QUESO::GslMatrix, // type for parameter matrices
         QUESO::GslVector, // type for qoi vectors
         QUESO::GslMatrix  // type for qoi matrices
        >(*env);

  //************************************************
  // Finalize environment
  //************************************************
  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
}
