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

#ifndef UQ_TEUCHOS_VECTOR_H
#define UQ_TEUCHOS_VECTOR_H

#ifdef QUESO_HAS_TRILINOS
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_LAPACK.hpp>
#endif

#include <queso/Defines.h>
#include <queso/Vector.h>
#include <vector>

#ifdef QUESO_HAS_TRILINOS

namespace QUESO {

/*! \file TeuchosVector.h
    \brief Vector class using Trilinos Teuchos
*/

/*! \class TeuchosVector
    \brief Class for vector operations using Teuchos (Trilinos).

    This class creates and provides basic support for vectors of templated
    type as a specialization of Vector using Teuchos vectors (from Trilinos), which are defined
    by an encapsulated Teuchos::SerialDenseVector structure.
*/


class TeuchosVector : public Vector
{
public:

  //! @name Constructor/Destructor methods.
  //@{

  //! Shaped constructor: creates an empty vector of size given by Map& map.
  TeuchosVector(const BaseEnvironment& env, const Map& map);

 //! Shaped constructor: creates an vector of values \c value and of size given by Map& map.
  TeuchosVector(const BaseEnvironment& env, const Map& map, double value);

 //! Shaped constructor: creates an vector of size given by Map& map and of values given by an average involving \c d1, \c d2 and the vector size.
  TeuchosVector(const BaseEnvironment& env, double d1, double d2, const Map& map);

  //! Shaped constructor: creates an vector of size given by vector \c v and of values given by an average involving \c d1, \c d2 and \c \this vector size.
  TeuchosVector(const TeuchosVector&         v, double d1, double d2);

  //! Copy constructor.
  TeuchosVector(const TeuchosVector&         y);

  //! Destructor
  ~TeuchosVector();
   //@}


  //! @name Set methods.
  //@{
  //! Set all values in the vector to a constant value.
  TeuchosVector& operator= (double a);

  //! Copies values from one vector to another.
  TeuchosVector& operator= (const TeuchosVector& rhs);

  //! Stores in \c this vector the coordinate-wise multiplication of \c this and a.
  TeuchosVector& operator*=(double a);

  //! Stores in \c this vector the coordinate-wise division of \c this by a.
  TeuchosVector& operator/=(double a);

  //! Stores in \c this vector the coordinate-wise multiplication of \c this with rhs.
  TeuchosVector& operator*=(const TeuchosVector& rhs);

  //! Stores in \c this vector the coordinate-wise division of \c this by rhs.
  TeuchosVector& operator/=(const TeuchosVector& rhs);

   //! Stores in \c this vector the coordinate-wise addition of \c this and rhs.
  TeuchosVector& operator+=(const TeuchosVector& rhs);

   //! Stores in \c this vector the coordinate-wise subtraction of \c this and \c rhs.
  TeuchosVector& operator-=(const TeuchosVector& rhs);
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

  //! @name Norm methods
  //@{
  //! Returns the norm of the vector, as the square root of 2-norm of this vector.
  double       norm2Sq          () const;

  //! Returns the 2-norm (Euclidean norm) of the vector.
  double       norm2            () const;

  //! Returns the 1-norm of the vector.
  double       norm1            () const;

    //! Returns the infinity-norm (maximum norm) of the vector.
  double       normInf          () const;


  //@}


  //! @name Comparison methods.
  //@{
  //! This function returns true if at least one component of \c this is smaller than the respective component of rhs.
  bool         atLeastOneComponentSmallerThan       (const TeuchosVector& rhs) const;

  //! This function returns true if at least one component of \c this is bigger than the respective component of rhs.
  bool         atLeastOneComponentBiggerThan        (const TeuchosVector& rhs) const;

  //! This function returns true if at least one component of \c this is smaller than or equal to the respective component of rhs.
  bool         atLeastOneComponentSmallerOrEqualThan(const TeuchosVector& rhs) const;

  //! This function returns true if at least one component of \c this is bigger than or equal to the respective component of rhs.
  bool         atLeastOneComponentBiggerOrEqualThan (const TeuchosVector& rhs) const;
  //@}

    //! @name Set methods.
  //@{
  //! Component-wise sets all values to \c this with value.
  void	       cwSet(double value);

  //! Component-wise sets all values of \c this with vector \c vec, starting at position \c initialPos.
  void         cwSet(unsigned int initialPos, const TeuchosVector& vec);

  //! Component-wise extracts all values of \c this with vector \c vec, starting at position \c initialPos.
  void         cwExtract(unsigned int initialPos, TeuchosVector& vec) const;

  //! This function inverts component-wise the element values of \c this.
  void         cwInvert();

  //! Component-wise sets the square-root of \c this.
  void         cwSqrt();

  //! This function concatenates vectors \c v1 and \c v2 into \c this vector.
  void         cwSetConcatenated(const TeuchosVector& v1, const TeuchosVector& v2);

  //! This function sets component-wise Gaussian random variates, with mean \c mean and standard deviation \c stdDev.
  void	       cwSetGaussian(double mean, double stdDev);

  //! This function sets component-wise Gaussian random variates, with vectors for mean and standard deviation.
  void	       cwSetGaussian(const TeuchosVector& meanVec,       const TeuchosVector& stdDevVec);

  //! This function sets component-wise a number uniformly distributed in the range of elements of [lowerBoundVec,upperBoundVec].
  void         cwSetUniform (const TeuchosVector& lowerBoundVec, const TeuchosVector& upperBoundVec);

  //! This function sets component-wise random variates from the Beta distribution, with vector parameters alpha and beta.
  void         cwSetBeta        (const TeuchosVector& alpha,   const TeuchosVector& beta     );

  //! This function sets component-wise random variates from the Inverse Gamma distribution, with parameters given by vectors  \c a and \c b.
  void         cwSetGamma       (const TeuchosVector& a,       const TeuchosVector& b        );

  //! This function sets component-wise random variates from the Inverse Gamma distribution, with parameters given by vectors  \c a and \c b.
  void         cwSetInverseGamma(const TeuchosVector& a,   const TeuchosVector& b     );
  //@}

   //! @name Miscellaneous methods.
  //@{

   //! This function returns absolute value of all elements in \c this.
  TeuchosVector abs() const;

  //! Copies the values of \c this vector (a TeuchosVector) to a std::vector structure.
  /*! With this method, the std::vector copy of the TeuchosVector may be sorted using std::sort.*/
  void		copy_to_std_vector(std::vector<double>& vec); // output

  //! Copies the values of std::vector structure to \c this vector (a TeuchosVector).
  /*! With this method, after a std::vector is sorted, it may be copied to a TeuchosVector.*/
  void		copy_from_std_vector(const std::vector<double> vec);

  //! This function sorts the elements of the vector \c this in ascending numerical order.
  /*! It copies \c this vector to a std::vector and uses the std::sort functionality to sort it. */
  void         sort             ();

  //! Returns the sum of the components of the vector.
  double       sumOfComponents  () const;

  //! Broadcasts a message from the process with \c srcRank root to all other processes of the group.
  void         mpiBcast         (int srcRank, const MpiComm& bcastComm);

  //! Combines values from all processes and distributes the result back to all processes.
  void         mpiAllReduce     (RawType_MPI_Op mpiOperation, const MpiComm& opComm, TeuchosVector& resultVec) const;

  //!  Gathers values from a group of processes and returns all quantiles.
  void         mpiAllQuantile   (double probability, const MpiComm& opComm, TeuchosVector& resultVec) const;

  //!  Reproduces MATLAB linear inter/extra-polation.
  /*! yiVec.matlabLinearInterpExtrap(xVec,yVec,xiVec) interpolates/extrapolates to find yiVec, the values of the
   * underlying function yVec at the points in the vector xiVec. Thus, yVec and xVec must have the same size;
   * and yiVec and xiVec must have the same size.*/
  void         matlabLinearInterpExtrap(const TeuchosVector& xVec, const TeuchosVector& yVec, const TeuchosVector& xiVec);

  //TODO
  void         matlabDiff       (unsigned int firstPositionToStoreDiff, double valueForRemainderPosition, TeuchosVector& outputVec) const;

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





private:
  //! Default Constructor
  /*! Creates an empty vector of no length.*/
  TeuchosVector();

  //! Teuchos vector.
  Teuchos::SerialDenseVector<int,double> m_vec;

  //! This function copies the elements of the vector \c src into \c this.
  void         copy             (const TeuchosVector& src);

};

//TeuchosVector copy		  (int , 						    double []);
TeuchosVector operator/    (double a,            		    const TeuchosVector& x  );
TeuchosVector operator/    (const TeuchosVector& x,   const TeuchosVector& y  );
TeuchosVector operator*    (double a,		    	        const TeuchosVector& x  );
TeuchosVector operator*    (const TeuchosVector& x,   const TeuchosVector& y  );
double               scalarProduct(const TeuchosVector& x,   const TeuchosVector& y  );
TeuchosVector operator+    (const TeuchosVector& x,   const TeuchosVector& y  );
TeuchosVector operator-    (const TeuchosVector& x,   const TeuchosVector& y  );
bool                 operator==   (const TeuchosVector& lhs, const TeuchosVector& rhs);
std::ostream&        operator<<   (std::ostream& os,                const TeuchosVector& obj);

}  // End namespace QUESO

#endif // ifdef QUESO_HAS_TRILINOS

#endif //UQ_TEUCHOS_VECTOR_H
