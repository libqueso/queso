//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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

#ifndef __UQ_GSL_VECTOR_H__
#define __UQ_GSL_VECTOR_H__

#include <uqVector.h>
#include <gsl/gsl_vector.h>

class uqGslVectorClass : public uqVectorClass
{
public:
  uqGslVectorClass();
  uqGslVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map);
  uqGslVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double value);
  uqGslVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const uqMapClass& map); // MATLAB linspace
  uqGslVectorClass(const uqGslVectorClass&         v, double d1, double d2);                        // MATLAB linspace
  uqGslVectorClass(const uqGslVectorClass&         y);
 ~uqGslVectorClass();

  uqGslVectorClass& operator= (const uqGslVectorClass& rhs);
  uqGslVectorClass& operator*=(double a);
  uqGslVectorClass& operator/=(double a);
  uqGslVectorClass& operator*=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator/=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator+=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator-=(const uqGslVectorClass& rhs);
            double& operator[](unsigned int i);
      const double& operator[](unsigned int i) const;

  unsigned int sizeLocal        () const;
  unsigned int sizeGlobal       () const;
  double       norm2Sq          () const;
  double       norm2            () const;
  double       norm1            () const;
  double       normInf          () const;
  double       sumOfComponents  () const;
  void         cwSet            (double value);
  void         cwSetGaussian    (const gsl_rng* rng, double mean, double stdDev);
  void         cwSetGaussian    (const gsl_rng* rng, const uqGslVectorClass& meanVec, const uqGslVectorClass& stdDevVec);
  void         cwSetUniform     (const gsl_rng* rng, const uqGslVectorClass& aVec,    const uqGslVectorClass& bVec     );
  void         cwSetBeta        (const gsl_rng* rng, const uqGslVectorClass& alpha,   const uqGslVectorClass& beta     );
  void         cwSetGamma       (const gsl_rng* rng, const uqGslVectorClass& a,       const uqGslVectorClass& b        );
  void         cwSetInverseGamma(const gsl_rng* rng, const uqGslVectorClass& alpha,   const uqGslVectorClass& beta     );
  void         cwSetConcatenated(const uqGslVectorClass& v1, const uqGslVectorClass& v2);
  void         cwSetConcatenated(const std::vector<const uqGslVectorClass*>& vecs);
  void         cwSet            (unsigned int initialPos, const uqGslVectorClass& vec);
  void         cwExtract        (unsigned int initialPos, uqGslVectorClass& vec) const;
  void         cwInvert         ();
  void         cwSqrt           ();
  void         matlabDiff       (unsigned int firstPositionToStoreDiff, double valueForRemainderPosition, uqGslVectorClass& outputVec) const;
  void         sort             ();
  void         mpiBcast         (int srcRank, const uqMpiCommClass& bcastComm);
  void         mpiAllReduce     (uqRawType_MPI_Op mpiOperation, const uqMpiCommClass& opComm, uqGslVectorClass& resultVec) const;
  void         mpiAllQuantile   (double probability, const uqMpiCommClass& opComm, uqGslVectorClass& resultVec) const;
  void         print            (std::ostream& os) const;
  void         subWriteContents (const std::string&            varNamePrefix,
                                 const std::string&            fileName,
                                 const std::string&            fileType,
                                 const std::set<unsigned int>& allowedSubEnvIds) const;
  void         subReadContents  (const std::string&            fileName,
                                 const std::string&            fileType,
                                 const std::set<unsigned int>& allowedSubEnvIds);

  bool         atLeastOneComponentSmallerOrEqualThan(const uqGslVectorClass& rhs) const;
  bool         atLeastOneComponentBiggerOrEqualThan (const uqGslVectorClass& rhs) const;

  // Necessary for uqGslMatrixClass::invertMultiply() and uqGslMatrixClass::setRow/Column
  gsl_vector*  data                          () const; 

  double       getMaxValue      () const;
  double       getMinValue      () const;
  int          getMaxValueIndex () const;
  int          getMinValueIndex () const;
  void         getMaxValueAndIndex( double& value, int& index );
  void         getMinValueAndIndex( double& value, int& index );

  uqGslVectorClass abs() const;

private:

  void         copy             (const uqGslVectorClass& src);

  gsl_vector* m_vec;
};

uqGslVectorClass operator/    (      double a,              const uqGslVectorClass& x  );
uqGslVectorClass operator/    (const uqGslVectorClass& x,   const uqGslVectorClass& y  );
uqGslVectorClass operator*    (      double a,              const uqGslVectorClass& x  );
uqGslVectorClass operator*    (const uqGslVectorClass& x,   const uqGslVectorClass& y  );
double           scalarProduct(const uqGslVectorClass& x,   const uqGslVectorClass& y  );
uqGslVectorClass operator+    (const uqGslVectorClass& x,   const uqGslVectorClass& y  );
uqGslVectorClass operator-    (const uqGslVectorClass& x,   const uqGslVectorClass& y  );
bool             operator==   (const uqGslVectorClass& lhs, const uqGslVectorClass& rhs);
std::ostream&    operator<<   (std::ostream& os,            const uqGslVectorClass& obj);

#endif // __UQ_GSL_VECTOR_H__
