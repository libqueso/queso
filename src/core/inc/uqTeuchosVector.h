//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012 The PECOS Development Team
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

#ifndef __UQ_TEUCHOS_VECTOR_H__
#define __UQ_TEUCHOS_VECTOR_H__

#ifdef QUESO_HAS_TRILINOS
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_LAPACK.hpp>
#endif
#include <vector>
#include <uqVector.h>
#include "boost/random.hpp"
#include <boost/math/distributions.hpp>  //for beta distributions

#ifdef QUESO_HAS_TRILINOS

class uqTeuchosVectorClass : public uqVectorClass
{
public:
  uqTeuchosVectorClass();
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map);
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double value);
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const uqMapClass& map);
  uqTeuchosVectorClass(const uqTeuchosVectorClass&         v, double d1, double d2);                   
  uqTeuchosVectorClass(const uqTeuchosVectorClass&         y);
  
  ~uqTeuchosVectorClass();

  uqTeuchosVectorClass abs() const;

  uqTeuchosVectorClass& operator= (double a);
  uqTeuchosVectorClass& operator= (const uqTeuchosVectorClass& rhs);
  uqTeuchosVectorClass& operator*=(double a);
  uqTeuchosVectorClass& operator/=(double a);
  uqTeuchosVectorClass& operator*=(const uqTeuchosVectorClass& rhs);
  uqTeuchosVectorClass& operator/=(const uqTeuchosVectorClass& rhs);
  uqTeuchosVectorClass& operator+=(const uqTeuchosVectorClass& rhs);
  uqTeuchosVectorClass& operator-=(const uqTeuchosVectorClass& rhs);
  

  double& operator[](unsigned int i);
  const double& operator[](unsigned int i) const;
  
  unsigned int sizeLocal        () const;
  unsigned int sizeGlobal       () const;
  double*      values           () ;
  double       norm2Sq          () const;
  double       norm2            () const;
  double       norm1            () const;
  double       normInf          () const;
  double       sumOfComponents  () const;
  
  double       getMaxValue      () const;
  double       getMinValue      () const;
  int          getMaxValueIndex () const;
  int          getMinValueIndex () const;
  void         getMaxValueAndIndex( double& value, int& index );
  void         getMinValueAndIndex( double& value, int& index );
  
  bool         atLeastOneComponentSmallerThan       (const uqTeuchosVectorClass& rhs) const;
  bool         atLeastOneComponentBiggerThan        (const uqTeuchosVectorClass& rhs) const;
  bool         atLeastOneComponentSmallerOrEqualThan(const uqTeuchosVectorClass& rhs) const;
  bool         atLeastOneComponentBiggerOrEqualThan (const uqTeuchosVectorClass& rhs) const;
  
  void	       cwSet(double value);
  void         cwSet(unsigned int initialPos, const uqTeuchosVectorClass& vec);
  void         cwExtract(unsigned int initialPos, uqTeuchosVectorClass& vec) const; 
  void         cwInvert();
  void         cwSqrt();
  void	       cwSetGaussian(double mean, double stdDev);
  void	       cwSetGaussian(const uqTeuchosVectorClass& meanVec,       const uqTeuchosVectorClass& stdDevVec);
  void         cwSetUniform (const uqTeuchosVectorClass& lowerBoundVec, const uqTeuchosVectorClass& upperBoundVec);
  void         cwSetUniform2(const uqTeuchosVectorClass& lowerBoundVec, const uqTeuchosVectorClass& upperBoundVec);//dummy
   
  void         cwSetConcatenated(const uqTeuchosVectorClass& v1, const uqTeuchosVectorClass& v2);
  void         cwSetBeta        (const uqTeuchosVectorClass& alpha,   const uqTeuchosVectorClass& beta     );
  void         cwSetInverseGamma(const uqTeuchosVectorClass& alpha,   const uqTeuchosVectorClass& beta     );
  
  void         cwSetGaussian2   (double mean, double stdDev);
  void         cwSetGaussian2   (const uqTeuchosVectorClass& meanVec, const uqTeuchosVectorClass& stdDevVec);
  void         cwSetGamma       (const uqTeuchosVectorClass& a,       const uqTeuchosVectorClass& b        );
  void         sort             ();
  void         print            (std::ostream& os) const;
  void         subReadContents  (const std::string&            fileName,
								 const std::string&            fileType,
								 const std::set<unsigned int>& allowedSubEnvIds);
  void         subWriteContents (const std::string&            varNamePrefix,
                                 const std::string&            fileName,
                                 const std::string&            fileType,
                                 const std::set<unsigned int>& allowedSubEnvIds) const;
				 
  double		GetRandomDoubleUsingNormalDistribution(int seed, double mean,double sigma);
  double	 	GetRandomDoubleUsingUniformZeroOneDistribution(int seed);
  void 			copy_to_std_vector(std::vector<double>& vec); // output
  void 			copy_from_std_vector(const std::vector<double> vec); //(input, output)


private:
  Teuchos::SerialDenseVector<int,double> m_vec;
  void         copy             (const uqTeuchosVectorClass& src);

};

uqTeuchosVectorClass copy		  (int , 						    double []);
uqTeuchosVectorClass operator/    (double a,            		    const uqTeuchosVectorClass& x  );
uqTeuchosVectorClass operator/    (const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
uqTeuchosVectorClass operator*    (double a,		    	        const uqTeuchosVectorClass& x  );
uqTeuchosVectorClass operator*    (const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
double               scalarProduct(const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
uqTeuchosVectorClass operator+    (const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
uqTeuchosVectorClass operator-    (const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
bool                 operator==   (const uqTeuchosVectorClass& lhs, const uqTeuchosVectorClass& rhs);
std::ostream&        operator<<   (std::ostream& os,                const uqTeuchosVectorClass& obj);

#endif // ifdef QUESO_HAS_TRILINOS

#endif //__UQ_TEUCHOS_VECTOR_H__


