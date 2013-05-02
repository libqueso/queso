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

#ifndef __UQ_TEUCHOS_MATRIX_H__
#define __UQ_TEUCHOS_MATRIX_H__

#ifdef QUESO_HAS_TRILINOS
#include <uqTeuchosVector.h>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_LAPACK.hpp>
#endif

#include <uqMatrix.h>

#ifdef QUESO_HAS_TRILINOS

/*! \file uqTeuchosMatrix.h
    \brief Vector class using Trilinos Teuchos
*/

/*! \class uqTeuchosMatrixClass
    \brief Class for matrices operations using Teuchos (Trilinos).
    
    This class creates and provides basic support for matrices of templated 
    type as a specialization of uqMatrixClass using Teuchos matrices (from Trilinos), which are defined 
    by an encapsulated Teuchos::SerialDenseMatrix structure.
*/

class uqTeuchosMatrixClass : public uqMatrixClass
{
public:
   //! @name Constructor/Destructor methods
  //@{ 

  //! Default Constructor
  /*! Creates an empty matrix vector of no dimension. It should not be used by user.*/
  uqTeuchosMatrixClass();
  
   //! Shaped Constructor: creates a shaped matrix with \c numCols columns.  
  uqTeuchosMatrixClass(const uqBaseEnvironmentClass& env,
                       const uqMapClass&             map,
                       unsigned int                  numCols);
  
    //! Shaped Constructor: creates a square matrix with size \c map.NumGlobalElements() and diagonal values all equal to \c diagValue.  
  uqTeuchosMatrixClass(const uqBaseEnvironmentClass& env,
                       const uqMapClass&             map,
     	               double                        diagValue); 
  
    //! Shaped Constructor: creates a square matrix with size \c v.sizeLocal() and diagonal values all equal to \c diagValue.  
  uqTeuchosMatrixClass(const uqTeuchosVectorClass&   v,
                   	   double                    diagValue); 
  
  //! Shaped Constructor: creates a square matrix with size \c v.sizeLocal().  
  /*! The diagonal values of this matrix are the elements in vector \c v. */  
  uqTeuchosMatrixClass(const uqTeuchosVectorClass&   v);         
  
  //! Shaped Constructor: \c this matrix is a copy of matrix \c B.
  uqTeuchosMatrixClass(const uqTeuchosMatrixClass&   B);
  
    //! Destructor
 ~uqTeuchosMatrixClass();
//@}

  //! @name Set methods
  //@{ 
  //! Copies values from matrix \c rhs to \c this. 
  uqTeuchosMatrixClass& operator= (const uqTeuchosMatrixClass& rhs);
  
  //! Stores in \c this the coordinate-wise multiplication of \c this and \c a.
  uqTeuchosMatrixClass& operator*=(double a);
  
  //! Stores in \c this the coordinate-wise division of \c this by \c a.
  uqTeuchosMatrixClass& operator/=(double a);
  
  //! Stores in \c this the coordinate-wise addition of \c this and \c rhs.
  uqTeuchosMatrixClass& operator+=(const uqTeuchosMatrixClass& rhs);
  //! Stores in \c this the coordinate-wise subtraction of \c this by \c rhs.
  uqTeuchosMatrixClass& operator-=(const uqTeuchosMatrixClass& rhs);
  //@}
  
    //! @name Accessor methods
  //@{
  //! Element access method (non-const).
  double& operator()(unsigned int i, unsigned int j);
  
  //! Element access method (const).  
  const double& operator()(unsigned int i, unsigned int j) const;
  //@}

  //! @name Attribute methods
  //@{ 
  //! Returns the local row dimension of \c this matrix.
  unsigned int      numRowsLocal           () const;
  
  //! Returns the global row dimension of \c this matrix.
  unsigned int      numRowsGlobal          () const;
  
  //! Returns the column dimension of \c this matrix.
  unsigned int      numCols                () const;
  
  //! Returns a pointer to the first element of \c this matrix.
  double*           values                 () ;//added by Kemelli on 12/04/12
  
  //! Returns the stride between the columns of this matrix in memory. 
  int 	            stride                 () ;//added by Kemelli on 12/04/12
  
  //! Returns the maximum element value of the matrix.
  double            max                       () const;
  
  //! This function returns the number of singular values of \c this matrix (rank). 
  /*! The rank function provides an estimate of the number of linearly independent rows or columns of a full matrix. */
  unsigned int      rank                      (double absoluteZeroThreshold, double relativeZeroThreshold) const;
  
  //! This function calculates the transpose of \c this matrix  (square).
  uqTeuchosMatrixClass  transpose                 () const;
  
  //! This function calculates the inverse of \c this matrix  (square).
  uqTeuchosMatrixClass  inverse                   () const;
  
  //! Calculates the determinant of \c this matrix.
  double                determinant               () const;
  
    //! Calculates the ln(determinant) of \c this matrix.
  double                lnDeterminant             () const;
  //@}
  
  //! @name Norm methods
  //@{ 
  //! Returns the Frobenius norm of \c this matrix.
  double            normFrob                  () const;
    
  //! Returns the Frobenius norm of \c this matrix.
  double            normMax                   () const;
  //@}
  
  //! @name Mathematical methods
  //@{ 
  //! Computes Cholesky factorization of  \c this, a real symmetric positive definite matrix. 
  /*! In case \this fails to be symmetric and positive definite, an error will be returned. */
  int               chol                   () ; 
  
  //! Checks for the dimension of \c this matrix, \c matU, \c VecS and \c matVt, and calls the protected routine \c internalSvd to compute the singular values of \c this. 
  int               svd                    (uqTeuchosMatrixClass& matU, uqTeuchosVectorClass& vecS, uqTeuchosMatrixClass& matVt) const;
  
  //! This function calls private member uqTeuchosMatrixClass::internalSvd() to set a M-by-N orthogonal matrix U of the singular value decomposition (svd) of a general rectangular M-by-N matrix A. 
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  const uqTeuchosMatrixClass& svdMatU      () const;
  
  //! This function calls private member  uqTeuchosMatrixClass::internalSvd() to set a N-by-N orthogonal square matrix V of the singular value decomposition (svd) of a general rectangular M-by-N matrix A. 
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  const uqTeuchosMatrixClass& svdMatV      () const;
  
  //! This function solves the system A x = b using the singular value decomposition (U, S, V) of A which must have been computed previously with uqTeuchosMatrixClass::svd (x=solVec, b=rhsVec). 
 /*! An orthogonal matrix A has a norm-preserving property, i.e. for any vector v, ||Av|| = ||v||. Then: 
 min(||Ax − b||^2) = min(||Ax − b||) = min(||UDVT x − b||) = min(||DV x − U b||).
 Substituting y = VT x and b' = UT b gives us Dy = b' with D being a diagonal matrix. 
 Or, y = inv(D)*UT*b and we only have to solve the linear system: VT x = y. */
  int               svdSolve               (const uqTeuchosVectorClass& rhsVec, uqTeuchosVectorClass& solVec) const;
  
  //! This function solves the system A x = b using the singular value decomposition (U, S, V) of A which must have been computed previously with uqTeuchosMatrixClass::svd (x=solMat, b=rhsMat).
  /*! Note that solMat is a matrix. Thus, to solve the system of equations '\this * solMath = rhsMat', this method calls 
   * svdSolve(const uqTeuchosVectorClass& rhsVec, uqTeuchosVectorClass& solVec)
   * passing one column of solMat with its repective column of rhsMat, .*/
  int               svdSolve               (const uqTeuchosMatrixClass& rhsMat, uqTeuchosMatrixClass& solMat) const;

  
  
  //! This function multiplies \c this matrix by vector \c x and returns a vector.
  uqTeuchosVectorClass  multiply                  (const uqTeuchosVectorClass& x) const;
  
  //! This function calculates the inverse of \c this matrix, multiplies it with vector \c b and stores the result in vector \c x.
  /*! It checks for a previous LU decomposition of \c this matrix and does not recompute it
   if private attribute m_LU != NULL .*/
  void                  invertMultiply            (const uqTeuchosVectorClass& b, uqTeuchosVectorClass& x) const;
  
  //! This function calculates the inverse of \c this matrix and multiplies it with vector \c b. 
  /*! It calls void uqTeuchosMatrixClass::invertMultiply(const uqTeuchosVectorClass& b, uqTeuchosVectorClass& x) internally.*/
  uqTeuchosVectorClass  invertMultiply            (const uqTeuchosVectorClass& b) const;
  
  //! This function calculates the inverse of \c this matrix, multiplies it with matrix \c B and stores the result in matrix \c X.
  /*! It checks for a previous LU decomposition of \c this matrix and does not recompute it if private attribute \c m_LU != NULL .*/
  void                  invertMultiply            (const uqTeuchosMatrixClass& B, uqTeuchosMatrixClass& X) const;
  
  //! This function calculates the inverse of \c this matrix and multiplies it with matrix \c B.
  /*! It calls void uqTeuchosMatrixClass::invertMultiply(const uqTeuchosMatrixClass& B, uqTeuchosMatrixClass& X) const internally.*/
  uqTeuchosMatrixClass  invertMultiply            (const uqTeuchosMatrixClass& B) const;
  
  //! This function calculates the inverse of \c this matrix, multiplies it with vector \c b and stores the result in vector \c x.
  /*! It recalculates the LU decomposition of \c this matrix.*/
  void                  invertMultiplyForceLU     (const uqTeuchosVectorClass& b, uqTeuchosVectorClass& x) const;
  
  //! This function calculates the inverse of \c this matrix and multiplies it with vector \c b. 
  /*! It calls void uqTeuchosMatrixClass::InvertMultiplyForceLU(const uqTeuchosVectorClass& b, uqTeuchosVectorClass& x) const;(const uqTeuchosVectorClass& b,
uqTeuchosVectorClass& x) internally.*/
  uqTeuchosVectorClass  invertMultiplyForceLU     (const uqTeuchosVectorClass& b) const;
    
  //! This function computes the eigenvalues of a real symmetric matrix.
  void              eigen                     (uqTeuchosVectorClass& eigenValues, uqTeuchosMatrixClass* eigenVectors) const;			
  
  //! This function finds largest eigenvalue, namely \c eigenValue, of \c this matrix and its corresponding eigenvector, namely \c eigenVector.
  void              largestEigen              (double& eigenValue, uqTeuchosVectorClass& eigenVector) const;
  
  //! This function finds smallest eigenvalue, namely \c eigenValue, of \c this matrix and its corresponding eigenvector, namely \c eigenVector.
  void              smallestEigen             (double& eigenValue, uqTeuchosVectorClass& eigenVector) const;
 
  //@}
  
  //! @name Get/Set methods
  //@{ 
      
  //! Component-wise set all values to \c this with value.
  void                  cwSet                     (double value);
  
  //! Set the components of \c which positions are greater than (rowId,colId) with the value of mat(rowId,colId).
  void                  cwSet                     (unsigned int rowId, unsigned int colId, const uqTeuchosMatrixClass& mat);

  
  //! This function gets the column_num-th column of \c this matrix and stores it into vector \c column.
  void                  getColumn                 (const unsigned int column_num, uqTeuchosVectorClass& column) const;
  
  //! This function gets the column_num-th column of \c this matrix.
  uqTeuchosVectorClass  getColumn                 (const unsigned int column_num) const;
  
  //! This function copies vector \c column into the column_num-th column of \c this matrix. 
  void                  setColumn                 (const unsigned int column_num, const uqTeuchosVectorClass& column);
  
  //! This function gets the row_num-th column of \c this matrix and stores it into vector \c row.
  void                  getRow                    (const unsigned int row_num, uqTeuchosVectorClass& row) const;
  
  //! This function gets the row_num-th column of \c this matrix.
  uqTeuchosVectorClass  getRow                    (const unsigned int row_num) const;
  
  //! This function copies vector \c row into the row_num-th column of \c this matrix. 
  void                  setRow                    (const unsigned int row_num, const uqTeuchosVectorClass& row);
  
  //! This function sets all the entries above the main diagonal  (inclusive or not) of \c this matrix to zero.
  /*! If \c includeDiagonal = false, then only the entries above the main diagonal are set to zero; 
  if \c includeDiagonal = true, then the elements of the matrix diagonal are also set to zero.*/
  void                  zeroLower              (bool includeDiagonal = false);
  
   //! This function sets all the entries bellow the main diagonal (inclusive or not) of \c this matrix to zero.
  /*! If \c includeDiagonal = false, then only the entries bellow the main diagonal are set to zero; 
  if \c includeDiagonal = true, then the elements of the matrix diagonal are also set to zero.*/
  void                  zeroUpper              (bool includeDiagonal = false);
  
  //! This function sets to zero (filters) all entries of \c this matrix which are smaller than \c thresholdValue.
  /*! If \c thresholdValue < 0 then no values will be filtered.*/
  void               filterSmallValues         (double thresholdValue);
  
  //! This function sets to zero (filters) all entries of \c this matrix which are greater than \c thresholdValue.
  /*! If \c thresholdValue < 0 then no values will be filtered.*/	
  void               filterLargeValues         (double thresholdValue);
  
   //! This function stores the transpose of \c this matrix into \c this matrix.
  void               fillWithTranspose         (const uqTeuchosMatrixClass& mat); 
  
  //! This function fills \c this matrix diagonally with const block  matrices.
  void               fillWithBlocksDiagonally  (const std::vector<const uqTeuchosMatrixClass* >& matrices);
  
  //! This function fills \c this matrix diagonally with block matrices.
  void               fillWithBlocksDiagonally  (const std::vector<      uqTeuchosMatrixClass* >& matrices);
  
  //! This function fills \c this matrix horizontally with const block matrices.
  void               fillWithBlocksHorizontally(const std::vector<const uqTeuchosMatrixClass* >& matrices);
  
  //! This function fills \c this matrix horizontally with block matrices.
  void               fillWithBlocksHorizontally(const std::vector<      uqTeuchosMatrixClass* >& matrices);
  
  //! This function fills \c this matrix vertically with const block matrices.
  void               fillWithBlocksVertically  (const std::vector<const uqTeuchosMatrixClass* >& matrices);
  
  //! This function fills \c this matrix vertically with block matrices.
  void               fillWithBlocksVertically  (const std::vector<      uqTeuchosMatrixClass* >& matrices);
  
  //! This function calculates the tensor product of matrices \c mat1 and \c mat2 and stores it in \c this matrix.
  void               fillWithTensorProduct     (const uqTeuchosMatrixClass& mat1, const uqTeuchosMatrixClass& mat2);
  
  //! This function calculates the tensor product of matrix \c mat1 and  vector \c vec2 and stores it in \c this matrix.
  void               fillWithTensorProduct     (const uqTeuchosMatrixClass& mat1, const uqTeuchosVectorClass& vec2); 
  //@}
  
 
    //! @name Miscellaneous methods
  //@{

  void              mpiSum                 (const uqMpiCommClass& comm, uqTeuchosMatrixClass& M_global) const;
  
  void              matlabLinearInterpExtrap  (const uqTeuchosVectorClass& x1Vec, const uqTeuchosMatrixClass& y1Mat, const uqTeuchosVectorClass& x2Vec);
  //@}
 
   //! @name I/O methods
  //@{
  //! Print method. Defines the behavior of the ostream << operator inherited from the Object class.   
  void                  print                     (std::ostream& os) const;
  
  //! Read contents of subenvironment from file \c fileName.
  void                  subReadContents           (const std::string&            fileName,
                                                   const std::string&            fileType,
                                                   const std::set<unsigned int>& allowedSubEnvIds);
  
  //! Write contents of subenvironment in file \c fileName.
  void                  subWriteContents          (const std::string&            varNamePrefix,
                                                   const std::string&            fileName,
                                                   const std::string&            fileType,
                                                   const std::set<unsigned int>& allowedSubEnvIds) const;
//@}						     
  
private:
  
  //! In this function \c this matrix receives a copy of matrix \c src.  
  void              copy                      (const uqTeuchosMatrixClass& src);
  
  //! In this function resets the LU decomposition of \c this matrix, as well as deletes the private member pointers, if existing.	
  void              resetLU                   ();
  
  //! This function multiplies \c this matrix by vector \c x and stores the resulting vector in \c y.
  void              multiply                  (const uqTeuchosVectorClass& x, uqTeuchosVectorClass& y) const;
	
  //! This function factorizes the M-by-N matrix A into the singular value decomposition A = U S V^T for M >= N. On output the matrix A is replaced by U.	
  /*! This function uses Teuchos GESVD computes the singular value decomposition (SVD) of a real  M-by-N matrix A, optionally computing 
   * the left and/or right singular vectors. The SVD is written A = U * SIGMA * transpose(V), where SIGMA is 
   * an M-by-N matrix which is zero except for its min(m,n) diagonal elements, U is an M-by-M orthogonal 
   * matrix, and V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA are the singular values 
   * of A; they are real and non-negative, and are returned in descending order.  The first min(m,n) columns 
   * of U and V are the left and right singular vectors of A. Note that the routine returns V**T, not V. */
  int               internalSvd               () const;

  //! Teuchos matrix, also referred to as \c this matrix.
  Teuchos::SerialDenseMatrix<int,double> m_mat;
  
  //! TECHOS matrix for the LU decomposition of m_mat.	  
  mutable Teuchos::SerialDenseMatrix<int,double> m_LU; 
  
  //! Stores the inverse of \c this matrix.
  mutable uqTeuchosMatrixClass* m_inverse;
  
  //! Mapping for matrices involved in the singular value decomposition (svd) routine.
  mutable uqMapClass*       	m_svdColMap;
  
  //! m_svdUmat stores the M-by-N orthogonal matrix U after the singular value decomposition of a matrix.
  mutable uqTeuchosMatrixClass* m_svdUmat;
  
  //! m_svdSvec stores the diagonal of the N-by-N diagonal matrix of singular values S after the singular value decomposition of a matrix.
  mutable uqTeuchosVectorClass* m_svdSvec;
  
  //! m_svdVmat stores the N-by-N orthogonal square matrix V after the singular value decomposition of a matrix.
  mutable uqTeuchosMatrixClass* m_svdVmat;
  
  //! m_svdVmatT stores the transpose of N-by-N orthogonal square matrix V, namely V^T, after the singular value decomposition of a matrix.
  mutable uqTeuchosMatrixClass* m_svdVTmat;
  
   //! The determinant of \c this matrix. 
  mutable double                m_determinant;
  
    //! The natural logarithm of the determinant of \c this matrix.
  mutable double                m_lnDeterminant;
  
  mutable uqTeuchosMatrixClass*  m_permutation;
  
  //! The pivoting vector of a LU decomposition.
  mutable int*              v_pivoting; 
  
  
  mutable int               m_signum;
  
  //! Indicates whether or not \c this matrix is singular.
  mutable bool              m_isSingular;
};

uqTeuchosMatrixClass operator*       (double a,                    const uqTeuchosMatrixClass& mat);
uqTeuchosVectorClass operator*       (const uqTeuchosMatrixClass& mat, const uqTeuchosVectorClass& vec);
uqTeuchosMatrixClass operator*       (const uqTeuchosMatrixClass& m1,  const uqTeuchosMatrixClass& m2 );
uqTeuchosMatrixClass operator+       (const uqTeuchosMatrixClass& m1,  const uqTeuchosMatrixClass& m2 );
uqTeuchosMatrixClass operator-       (const uqTeuchosMatrixClass& m1,  const uqTeuchosMatrixClass& m2 );
uqTeuchosMatrixClass matrixProduct   (const uqTeuchosVectorClass& v1,  const uqTeuchosVectorClass& v2 );
uqTeuchosMatrixClass leftDiagScaling (const uqTeuchosVectorClass& vec, const uqTeuchosMatrixClass& mat);
uqTeuchosMatrixClass rightDiagScaling(const uqTeuchosMatrixClass& mat, const uqTeuchosVectorClass& vec);
std::ostream&        operator<<      (std::ostream& os,            const uqTeuchosMatrixClass& obj);

#endif // ifdef QUESO_HAS_TRILINOS

#endif // __UQ_TEUCHOS_MATRIX_H__

