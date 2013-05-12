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

#include <uqBasicPdfsGsl.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>

// Default constructor ------------------------------
uqBasicPdfsGslClass::uqBasicPdfsGslClass()
  :
  uqBasicPdfsBaseClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_worldRank,
                      "uqBasicPdfsGslClass::constructor(), default",
                      "should not be used by user");
}

//! Constructor ---------------------------
uqBasicPdfsGslClass::uqBasicPdfsGslClass(int worldRank)
  :
  uqBasicPdfsBaseClass(worldRank)
{
}

// Destructor ---------------------------------------
uqBasicPdfsGslClass::~uqBasicPdfsGslClass()
{
}

// --------------------------------------------------
double
uqBasicPdfsGslClass::betaPdfActualValue(double x, double alpha, double beta) const
{
  return gsl_ran_beta_pdf(x,alpha,beta);
}

// --------------------------------------------------
double
uqBasicPdfsGslClass::gammaPdfActualValue(double x, double a, double b) const
{
  return gsl_ran_gamma_pdf(x,a,b);
}
