//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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
// $Id:$
//
//--------------------------------------------------------------------------

#include <example_likelihood.h>

double likelihoodRoutine(
  const uqGslVectorClass& paramValues,
  const uqGslVectorClass* paramDirection,
  const void*             functionDataPtr,
  uqGslVectorClass*       gradVector,
  uqGslMatrixClass*       hessianMatrix,
  uqGslVectorClass*       hessianEffect)
{
  //const uqGslVectorClass& meanVector =
  //  *((likelihoodRoutine_DataType *) functionDataPtr)->meanVector;
  //const uqGslMatrixClass& covMatrix  =
  //  *((likelihoodRoutine_DataType *) functionDataPtr)->covMatrix;

  //uqGslVectorClass diffVec(paramValues - meanVector);

  //return scalarProduct(diffVec, covMatrix.invertMultiply(diffVec));

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

#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
  double returnValue = -.5*resultValue;
#else
  double returnValue = resultValue; 
#endif

  if (paramValues.env().exceptionalCircunstance()) {
    if ((paramValues.env().subDisplayFile()       ) &&
        (paramValues.env().displayVerbosity() >= 0)) { // detailed output debug
      *paramValues.env().subDisplayFile() << "Leaving likelihood function"
                                          << ": paramValues = " << paramValues
                                          << ", returnValue = " << returnValue
                                          << std::endl;
    }
  }

  return returnValue;
}
