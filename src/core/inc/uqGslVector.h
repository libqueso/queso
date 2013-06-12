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
// Doxygen comments added by kemelli @ violeta Thu Nov  1 16:25:13 PDT 2012
//--------------------------------------------------------------------------

#ifndef __UQ_GSL_VECTOR_H__
#define __UQ_GSL_VECTOR_H__

/*! \file uqGslVector.h
    \brief Vector class using GSL
*/

#include <uqDefines.h>
#include <uqVector.h>
#include <gsl/gsl_vector.h>

/*! \class uqGslVectorClass

    \brief Class for vector operations using GSL library.
    
    This class creates and provides basic support for vectors of templated 
    type as a specialization of uqVectorClass using GSL vectors, which are defined 
    by an encapsulated gsl_vector structure.
*/

class uqGslVectorClass : public uqVectorClass
{
public:
  
  //! @name Constructor/Destructor methods.
  //@{ 

  //! Default Constructor
  /*! Creates an empty vector of no length.*/
  uqGslVectorClass();
  
  uqGslVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map);
  uqGslVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double value);
  uqGslVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const uqMapClass& map); // MATLAB linspace

  //! Construct a vector, with length the same as \c v, with evenly spaced numbers from \c start to \c end, inclusive
  uqGslVectorClass(const uqGslVectorClass&         v, double start, double end);
  uqGslVectorClass(const uqGslVectorClass&         y);
  
  //! Destructor
  ~uqGslVectorClass();
  //@}
 
  //! @name Set methods.
  //@{ 
  //! 	Copies values from vector rhs to \c this. 
  uqGslVectorClass& operator= (const uqGslVectorClass& rhs);
  
  //! Stores in \c this the coordinate-wise multiplication of \c this and a.
  uqGslVectorClass& operator*=(double a);
  
  //! Stores in \c this the coordinate-wise division of \c this by a.
  uqGslVectorClass& operator/=(double a);
  
  //! Stores in \c this the coordinate-wise multiplication of \c this with rhs.
  uqGslVectorClass& operator*=(const uqGslVectorClass& rhs);
  
  //! Stores in \c this the coordinate-wise division of \c this by rhs.
  uqGslVectorClass& operator/=(const uqGslVectorClass& rhs);
  
   //! Stores in \c this the coordinate-wise addition of \c this and rhs.
  uqGslVectorClass& operator+=(const uqGslVectorClass& rhs);
  
  //! Stores in \c this the coordinate-wise subtraction of \c this by rhs.
  uqGslVectorClass& operator-=(const uqGslVectorClass& rhs);
  //@}
  
  //! @name Accessor methods.
  //@{
  //! Element access method (non-const).
  /*! Returns the i-th element if x[i] is specified, the expression x(i) will return the same element.*/  
            double& operator[](unsigned int i);
	    
  //! Element access method (const).
  /*! Returns the i-th element if x[i] is specified, the expression x(i) will return the same element.*/ 
	    
      const double& operator[](unsigned int i) const;
  //@}

  //! @name Attribute methods.
  //@{ 
  //! Returns the length of this vector.
  unsigned int sizeLocal        () const;

  //! Returns the global length of this vector.
  unsigned int sizeGlobal       () const;
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
  
  //! @name Set methods.
  //@{ 
  //! Component-wise sets all values to \c this with value.
  void         cwSet            (double value);
  
  //! This function sets component-wise Gaussian random variates, with mean mean and standard deviation stdDev.
  void         cwSetGaussian    (double mean, double stdDev);
  
  //! This function sets component-wise Gaussian random variates, with vectors for mean and standard deviation.
  void         cwSetGaussian    (const uqGslVectorClass& meanVec, const uqGslVectorClass& stdDevVec);
  
  //! This function sets component-wise a number uniformly distributed in the range of elements of [aVec,bVec].
  void         cwSetUniform     (const uqGslVectorClass& aVec, const uqGslVectorClass& bVec);
  
  //! This function returns a random variate from the beta distribution, with vector parameters alpha and beta.
  void         cwSetBeta        (const uqGslVectorClass& alpha, const uqGslVectorClass& beta);
  
  //! This function returns a random variate from the gamma distribution with vector parameters a and b. 
  void         cwSetGamma       (const uqGslVectorClass& a, const uqGslVectorClass& b);
  
  //! This function returns a random variate from the inverse gamma distribution with vector parameters alpha and beta. 
  void         cwSetInverseGamma(const uqGslVectorClass& alpha, const uqGslVectorClass& beta);
  
  //! This function concatenates uqGslVectorClass v1 and uqGslVectorClass v2 into  \c this.
  void         cwSetConcatenated(const uqGslVectorClass& v1, const uqGslVectorClass& v2);
  
  //! This function concatenates vectors \c v1 and \c v2 into \c this vector.
  void         cwSetConcatenated(const std::vector<const uqGslVectorClass* >& vecs);
  
  //! This function sets the vector vec into \c this starting at position initialPos.  
  void         cwSet            (unsigned int initialPos, const uqGslVectorClass& vec);

  //! This function sets the values of this starting at position initialPos ans saves them in vector vec.  
  void         cwExtract        (unsigned int initialPos, uqGslVectorClass& vec) const;
  
  //! This function inverts component-wise the element values of \c this.  
  void         cwInvert         ();
  
  //! This function returns component-wise the square-root of \c this.  
  void         cwSqrt           ();
  //@}

//! @name I/O methods.
  //@{ 
  //! Print method.  Defines the behavior of the std::ostream << operator inherited from the Object class.
  void         print            (std::ostream& os) const;
  //@}

  //! This function sorts the elements of the vector \c this in ascending numerical order. 
  void         sort             ();
    
  void         matlabDiff       (unsigned int firstPositionToStoreDiff, double valueForRemainderPosition, uqGslVectorClass& outputVec) const;
  void         matlabLinearInterpExtrap(const uqGslVectorClass& x1Vec, const uqGslVectorClass& y1Vec, const uqGslVectorClass& x2Vec);
  void         mpiBcast         (int srcRank, const uqMpiCommClass& bcastComm);
  void         mpiAllReduce     (uqRawType_MPI_Op mpiOperation, const uqMpiCommClass& opComm, uqGslVectorClass& resultVec) const;
  void         mpiAllQuantile   (double probability, const uqMpiCommClass& opComm, uqGslVectorClass& resultVec) const;
  void         subWriteContents (const std::string&            varNamePrefix,
                                 const std::string&            fileName,
                                 const std::string&            fileType,
                                 const std::set<unsigned int>& allowedSubEnvIds) const;
  void         subReadContents  (const std::string&            fileName,
                                 const std::string&            fileType,
                                 const std::set<unsigned int>& allowedSubEnvIds);

//! @name Comparison methods.
  //@{ 
  //! This function returns true if at least one component of \c this is smaller than the respective component of rhs.
  bool         atLeastOneComponentSmallerThan       (const uqGslVectorClass& rhs) const;
  
  //! This function returns true if at least one component of \c this is bigger than the respective component of rhs. 
  bool         atLeastOneComponentBiggerThan        (const uqGslVectorClass& rhs) const;
  
  //! This function returns true if at least one component of \c this is smaller than or equal to the respective component of rhs.
  bool         atLeastOneComponentSmallerOrEqualThan(const uqGslVectorClass& rhs) const;
  
  //! This function returns true if at least one component of \c this is bigger than or equal to the respective component of rhs. 
  bool         atLeastOneComponentBiggerOrEqualThan (const uqGslVectorClass& rhs) const;
  //@}

  // Necessary for uqGslMatrixClass::invertMultiply() and uqGslMatrixClass::setRow/Column
  gsl_vector*  data                          () const; 

//! @name Attribute methods.
  //@{ 
  //! Returns the maximum value in the vector \c this.
  double       getMaxValue      () const;
  
  //! Returns minimum value in the vector \c this.
  double       getMinValue      () const;
  
  //! This function returns the index of the maximum value in the vector \c this.
  int          getMaxValueIndex () const;
  
  //! This function returns the index of the minimum value in the vector \c this.
  int          getMinValueIndex () const;
  
  //! This function returns maximum value in the vector \c this and its the index.
  void         getMaxValueAndIndex( double& value, int& index );
  
  //! This function returns minimum value in the vector \c this and its the index.
  void         getMinValueAndIndex( double& value, int& index );

  //! This function returns absolute value of elements in \c this.
  uqGslVectorClass abs() const;
 //@}
private:

  //! This function copies the elements of the vector src into \c this.
  void         copy             (const uqGslVectorClass& src);

  //! GSL vector.
  gsl_vector* m_vec;
};

// Comments in this part of file don't appear in the doxygen docs.

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
