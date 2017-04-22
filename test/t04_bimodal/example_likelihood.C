//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
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

#include <example_likelihood.h>
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
  //const uqGslVector& meanVector =
  //  *((likelihoodRoutine_DataType *) functionDataPtr)->meanVector;
  //const uqGslMatrix& covMatrix  =
  //  *((likelihoodRoutine_DataType *) functionDataPtr)->covMatrix;

  //uqGslVector diffVec(paramValues - meanVector);

  //return scalarProduct(diffVec, covMatrix.invertMultiply(diffVec));

  likelihoodCounter++;

  if (paramDirection  ||
      functionDataPtr ||
      gradVector      ||
      hessianMatrix   ||
      hessianEffect) {}; // just to remove compiler warning

  double x = paramValues[0];

  double mean1  = 10.;
  double sigma1 = 1.;
  double y1     = (x-mean1)*(x-mean1)/(2.*sigma1*sigma1);
//double z1     = (1./sigma1/sqrt(2*M_PI))*exp(-y1);
  double ln_z1  = log(1./sigma1/M_PI/sqrt(2*M_PI)) - y1;

  double mean2  = 100.;
  double sigma2 = 5.;
  double y2     = (x-mean2)*(x-mean2)/(2.*sigma2*sigma2);
//double z2     = (1./sigma2/sqrt(2*M_PI))*exp(-y2);
  double ln_z2  = log(1./sigma2/M_PI/sqrt(2*M_PI)) - y2;


//double resultValue = log((z1+2.*z2)/3.);
  double ln_zMin     = std::min(log(1./3.)+ln_z1, log(2./3.)+ln_z2);
  double ln_zMax     = std::max(log(1./3.)+ln_z1, log(2./3.)+ln_z2);
  double resultValue = log(1. + exp(ln_zMin - ln_zMax)) + ln_zMax;

  if (resultValue == -INFINITY) {
#if 1
    std::cerr << "WARNING In likelihoodRoutine"
              << ", likelihoodCounter =" << likelihoodCounter
              << ", fullRank "           << paramValues.env().fullRank()
              << ", subEnvironment "     << paramValues.env().subId()
              << ", subRank "            << paramValues.env().subRank()
              << ", inter0Rank "         << paramValues.env().inter0Rank()
              << ": x = "                << x
              << ", ln_z1 = "            << ln_z1
              << ", ln_z2 = "            << ln_z2
              << ", resultValue = "      << resultValue
              << std::endl;
#endif
    resultValue = -520.;
  }

  double returnValue = resultValue;

  if (paramValues.env().exceptionalCircumstance()) {
    if ((paramValues.env().subDisplayFile()      ) &&
        (paramValues.env().displayVerbosity() > 0)) { // detailed output debug
      *paramValues.env().subDisplayFile() << "Leaving likelihood function"
                                          << ": paramValues = " << paramValues
                                          << ", returnValue = " << returnValue
                                          << std::endl;
    }
  }

  return returnValue;
}
