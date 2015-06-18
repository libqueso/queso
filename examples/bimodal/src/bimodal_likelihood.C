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

#include <bimodal_likelihood.h>
#include <cmath>

static unsigned int likelihoodCounter = 0;

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  likelihoodCounter++;

  if (paramDirection  ||
      functionDataPtr ||
      gradVector      ||
      hessianMatrix   ||
      hessianEffect) {}; // just to remove compiler warning

    double returnValue = 0.;
    double x = paramValues[0];
    double mean1  = 10.;
    double sigma1 = 1.;
    double y1 = (x-mean1)*(x-mean1)/(2.*sigma1*sigma1);
    double z1 = (1./sigma1/sqrt(2*M_PI))*exp(-y1);

    double mean2  = 100.;
    double sigma2 = 5.;
    double y2 = (x-mean2)*(x-mean2)/(2.*sigma2*sigma2);
    double z2 = (1./sigma2/sqrt(2*M_PI))*exp(-y2);

    double resultValue = -2*log((z1+2.*z2)/3.);

    if (resultValue == INFINITY) {
      //std::cerr << "WARNING In likelihoodRoutine"
      //          << ", fullRank "       << paramValues.env().fullRank()
      //          << ", subEnvironment " << paramValues.env().subId()
      //          << ", subRank "        << paramValues.env().subRank()
      //          << ", inter0Rank "     << paramValues.env().inter0Rank()
      //          << ": x = "            << x
      //          << ", z1 = "           << z1
      //          << ", z2 = "           << z2
      //          << ", resultValue = "  << resultValue
      //          << std::endl;
      resultValue = 1040.;
    }

    returnValue = -.5*resultValue;


  if (paramValues.env().exceptionalCircumstance()) {
    if ((paramValues.env().subDisplayFile()       ) &&
        (paramValues.env().displayVerbosity() > 0)) { // detailed output debug
      *paramValues.env().subDisplayFile() << "Leaving likelihood function"
                                          << ": paramValues = " << paramValues
                                          << ", returnValue = " << returnValue
                                          << std::endl;
    }
  }

  return returnValue;
}
