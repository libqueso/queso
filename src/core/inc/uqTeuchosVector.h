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
// $Id:$
//
//--------------------------------------------------------------------------

#ifndef __UQ_TEUCHOS_VECTOR_H__
#define __UQ_TEUCHOS_VECTOR_H__

/*! \file uqTeuchosVector.h
    \brief Vector class using Trilinos Teuchos
*/

/*! \class uqTeuchosVectorClass
    \brief This class creates and provides basic support for vectors of templated 
    type as a specialization of uqVectorClass using Teuchos vectors (from Trilinos), which are defined 
    by an encapsulated Teuchos::SerialDenseVector structure.
*/

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
  
    //! @name Constructor/Destructor methods.
  //@{ 

  //! Default Constructor
  /*! Creates an empty vector of no length.*/
  uqTeuchosVectorClass();
  
  /*! Creates an empty vector of size given by uqMapClass& map.*/
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map);
  
  /*! Creates an vector of values \c value and of size given by uqMapClass& map.*/
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double value);
  
  /*! Creates an vector of size given by uqMapClass& map and of values given by an average envolving \c d1, \c d2 and the vector size.*/
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const uqMapClass& map);
 
  /*! Creates an vector of size given by vector \c v and of values given by an average envolving \c d1, \c d2 and \c \this vector size.*/
  uqTeuchosVectorClass(const uqTeuchosVectorClass&         v, double d1, double d2);                   
  
  /*! Copy constructor.*/
  uqTeuchosVectorClass(const uqTeuchosVectorClass&         y);
  
  //! Destructor
  
  ~uqTeuchosVectorClass();
   //@}
  



  //! @name Set methods.
  //@{ 
  //! Set all values in the vector to a constant value.  
  uqTeuchosVectorClass& operator= (double a);
  
  //! Copies values from one vector to another. 
  uqTeuchosVectorClass& operator= (const uqTeuchosVectorClass& rhs);
  
  //! Stores in \c this the cordinate-wise multiplication of \c this and a.
  uqTeuchosVectorClass& operator*=(double a);
  
  //! Stores in \c this the cordinate-wise division of \c this by a.
  uqTeuchosVectorClass& operator/=(double a);
  
  //! Stores in \c this the cordinate-wise multipication of \c this with rhs.
  uqTeuchosVectorClass& operator*=(const uqTeuchosVectorClass& rhs);
  
  //! Stores in \c this the cordinate-wise division of \c this by rhs.
  uqTeuchosVectorClass& operator/=(const uqTeuchosVectorClass& rhs);
  
   //! Stores in \c this the cordinate-wise addition of \c this and rhs.
  uqTeuchosVectorClass& operator+=(const uqTeuchosVectorClass& rhs);
  
   //! Stores in \c this the cordinate-wise subtraction of \c this and rhs.
  uqTeuchosVectorClass& operator-=(const uqTeuchosVectorClass& rhs);
  //@}

    //! @name Accessor methods.
  //@{
  //! Element access method (non-const).
  double& operator[](unsigned int i);
  
  //! Element access method (const).
  const double& operator[](unsigned int i) const;
  //@}

  //! @name Attribute methods.
  //@{ 
  //! Returns the length of this vector
  unsigned int sizeLocal        () const;
  
  //! Returns the global length of this vector.
  unsigned int sizeGlobal       () const;
  
  //TODO
  double*      values           () ;
  
  //@}
  
  //! @name Mathematical methods.
  //@{ 
  //! Returns the norm of the vector, as the square root of 2-norm of this vector.
  double       norm2Sq          () const;
  
  //! Returns the 2-norm (Euclidean norm) of the vector.
  double       norm2            () const;
  
  //! Returns the 1-norm of the vector.
  double       norm1            () const;
  
    //! Returns the infinity-norm (maximum norm) of the vector.
  double       normInf          () const;
  
    //! Returns the sum of the components of the vector.
  double       sumOfComponents  () const;
  //@} 
  
  //! @name Attribute methods.
  //@{ 
  //! Returns the maximum value in the \c this vector.
  double       getMaxValue      () const;
  
  //! Returns the minimum value in the \c this vector.
  double       getMinValue      () const;
  
  //! This function returns the index of the maximum value in the vector \c this.
  int          getMaxValueIndex () const;
  
  //! This function returns the index of the minimum value in the vector \c this.
  int          getMinValueIndex () const;
  
  //! This function returns maximum value in the vector \c this and its the index. 
  void         getMaxValueAndIndex( double& value, int& index );
  
    //! This function returns minimum value in the vector \c this and its the index.
  void         getMinValueAndIndex( double& value, int& index );
  //@}
  
  //! @name Comparison methods.
  //@{ 
  //! This function returns true if at least one component of \c this is smaller than the respective component of rhs.
  bool         atLeastOneComponentSmallerThan       (const uqTeuchosVectorClass& rhs) const;
  
  //! This function returns true if at least one component of \c this is bigger than the respective component of rhs. 
  bool         atLeastOneComponentBiggerThan        (const uqTeuchosVectorClass& rhs) const;
  
  //! This function returns true if at least one component of \c this is smaller than or equal to the respective component of rhs.
  bool         atLeastOneComponentSmallerOrEqualThan(const uqTeuchosVectorClass& rhs) const;
  
  //! This function returns true if at least one component of \c this is bigger than or equal to the respective component of rhs. 
  bool         atLeastOneComponentBiggerOrEqualThan (const uqTeuchosVectorClass& rhs) const;
  //@}
  
    //! @name Set methods.
  //@{ 
  //! Component-wise set all values to \c this with value.
  void	       cwSet(double value);
  
  //! Component-wise set all values of \c this with vector \c vec, starting at position \c initialPos
  void         cwSet(unsigned int initialPos, const uqTeuchosVectorClass& vec);
  
  //! Component-wise extract all values of \c this with vector \c vec, starting at position \c initialPos
  void         cwExtract(unsigned int initialPos, uqTeuchosVectorClass& vec) const; 
  
  //! This function inverts component-wise the element values of \c this.  
  void         cwInvert();
  
  //! Compoments-wise set the square-root of \c this.  
  void         cwSqrt();
  
  //! This function concatenates vectors \c v1 and \c v2 into \c this vector.
  void         cwSetConcatenated(const uqTeuchosVectorClass& v1, const uqTeuchosVectorClass& v2);
  
  //! This function sets component-wise Gaussian random variates, with mean \c mean and standard deviation \c stdDev. 
  void	       cwSetGaussian(double mean, double stdDev);
  
  //! This function sets component-wise Gaussian random variates, with vectors for mean and standard deviation.
  void	       cwSetGaussian(const uqTeuchosVectorClass& meanVec,       const uqTeuchosVectorClass& stdDevVec);
  
  //! This function sets component-wise a number uniformly distributed in the range of elements of [lowerBoundVec,upperBoundVec].
  void         cwSetUniform (const uqTeuchosVectorClass& lowerBoundVec, const uqTeuchosVectorClass& upperBoundVec);
  
  //void         cwSetUniform2(const uqTeuchosVectorClass& lowerBoundVec, const uqTeuchosVectorClass& upperBoundVec);//dummy
  
    //! This function sets component-wise random variates from the Beta distribution, with vector parameters alpha and beta.
  void         cwSetBeta        (const uqTeuchosVectorClass& alpha,   const uqTeuchosVectorClass& beta     );
  
  //! This function sets component-wise random variates from the Inverse Gamma distribution, with parameters given by vectors  \c a and \c b.
  void         cwSetGamma       (const uqTeuchosVectorClass& a,       const uqTeuchosVectorClass& b        );
    
  //! This function sets component-wise random variates from the Inverse Gamma distribution, with parameters given by vectors  \c a and \c b.
  void         cwSetInverseGamma(const uqTeuchosVectorClass& a,   const uqTeuchosVectorClass& b     );
  
  //void         cwSetGaussian2   (double mean, double stdDev);
  //void         cwSetGaussian2   (const uqTeuchosVectorClass& meanVec, const uqTeuchosVectorClass& stdDevVec);
  //@}
  
   //! @name Miscellaneous methods.
  //@{ 
    //! This function sorts the elements of the vector \c this in ascending numerical order. 
  void         sort             ();
  
   //! This function returns absolute value of elements in \c this.
  uqTeuchosVectorClass abs() const;
  //@}
  
  //! @name I/O methods.
  //@{ 
  //! Print method.  Defines the behavior of the std::ostream << operator inherited from the Object class.
  void         print            (std::ostream& os) const;
  void         subReadContents  (const std::string&            fileName,
								 const std::string&            fileType,
								 const std::set<unsigned int>& allowedSubEnvIds);
  void         subWriteContents (const std::string&            varNamePrefix,
                                 const std::string&            fileName,
                                 const std::string&            fileType,
                                 const std::set<unsigned int>& allowedSubEnvIds) const;
  //@} 
  void         mpiBcast         (int srcRank, const uqMpiCommClass& bcastComm);
  void         mpiAllReduce     (uqRawType_MPI_Op mpiOperation, const uqMpiCommClass& opComm, uqTeuchosVectorClass& resultVec) const;
  void         mpiAllQuantile   (double probability, const uqMpiCommClass& opComm, uqTeuchosVectorClass& resultVec) const;
 
  void         matlabLinearInterpExtrap(const uqTeuchosVectorClass& x1Vec, const uqTeuchosVectorClass& y1Vec, const uqTeuchosVectorClass& x2Vec);
  void         matlabDiff       (unsigned int firstPositionToStoreDiff, double valueForRemainderPosition, uqTeuchosVectorClass& outputVec) const;

                                 
				 
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


