//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#include <queso/BasicPdfsGsl.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>
#include <math.h>

namespace QUESO {

// Default constructor ------------------------------
BasicPdfsGsl::BasicPdfsGsl()
  :
  BasicPdfsBase()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      "BasicPdfsGsl::constructor(), default",
                      "should not be used by user");
}

//! Constructor ---------------------------
BasicPdfsGsl::BasicPdfsGsl(int worldRank)
  :
  BasicPdfsBase(worldRank)
{
}

// Destructor ---------------------------------------
BasicPdfsGsl::~BasicPdfsGsl()
{
}

// --------------------------------------------------
double
BasicPdfsGsl::betaPdfActualValue(double x, double alpha, double beta) const
{
  double result = gsl_ran_beta_pdf(x,alpha,beta);
  if (isinf(result)) { // CSRI - 2013-aug-06, with Laura
    std::cerr << "In BasicPdfsGsl::betaPdfActualValue(): hitting inf"
              << ", x = "     << x
              << ", alpha = " << alpha
              << ", beta = "  << beta
              << std::endl;
    //result = gsl_ran_beta_pdf(x,alpha,beta);
    result = 0.;
  }
  return result; 
}

// --------------------------------------------------
double
BasicPdfsGsl::gammaPdfActualValue(double x, double a, double b) const
{
  return gsl_ran_gamma_pdf(x,a,b);
}

}  // End namespace QUESO
