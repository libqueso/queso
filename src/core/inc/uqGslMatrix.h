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

#ifndef __UQ_GSL_MATRIX_H__
#define __UQ_GSL_MATRIX_H__

/*! \file uqGslMatrix.h
    \brief QUESO matrix class using GSL.
*/

#include <uqMatrix.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <uqGslVector.h>

/*! \class uqGslMatrixClass
    \brief Class for matrices operations using GSL library.
    
    This class creates and provides basic support for matrices of templated 
    type as a specialization of uqMatrixClass using GSL matrices, which are defined 
    by an encapsulated gsl_matrix structure.
*/

class uqGslMatrixClass : public uqMatrixClass
{
public:
 //! @name Constructor/Destructor methods
  //@{ 

  //! Default Constructor
  /*! Creates an empty matrix vector of no dimension. It should not be used by user.*/
  uqGslMatrixClass();
  
  //! Shaped Constructor: creates a shaped matrix with \c numCols columns.  
  uqGslMatrixClass(const uqBaseEnvironmentClass& env,
                   const uqMapClass&             map,
                   unsigned int                  numCols);
  
  //! Shaped Constructor: creates a square matrix with size \c map.NumGlobalElements() and diagonal values all equal to \c diagValue.  
  uqGslMatrixClass(const uqBaseEnvironmentClass& env,
                   const uqMapClass&             map,
                   double                        diagValue); // MATLAB eye
  
  //! Shaped Constructor: creates a square matrix with size \c v.sizeLocal() and diagonal values all equal to \c diagValue.  
  uqGslMatrixClass(const uqGslVectorClass&       v,
                   double                        diagValue); // MATLAB eye
  
  //! Shaped Constructor: creates a square matrix with size \c v.sizeLocal().  
  /*! The diagonal values of this matrix are the elements in vector \c v. */  
  uqGslMatrixClass(const uqGslVectorClass&       v);         // MATLAB diag
  
  //! Shaped Constructor: creates a matrix with  \c B.numCols() columns and \c B.numRowsLocal() rows.  
  /*! \c This matrix is a copy of matrix \c B. */  
  uqGslMatrixClass(const uqGslMatrixClass&       B);
  
  //! Destructor
  ~uqGslMatrixClass();
  //@}


  //! @name Set methods
  //@{ 
  //! 	Copies values from matrix \c rhs to \c this. 
  uqGslMatrixClass& operator= (const uqGslMatrixClass& rhs);
  
  //! Stores in \c this the coordinate-wise multiplication of \c this and \c a.
  uqGslMatrixClass& operator*=(double a);
  
  //! Stores in \c this the coordinate-wise division of \c this by \c a.
  uqGslMatrixClass& operator/=(double a);
  
  //! Stores in \c this the coordinate-wise addition of \c this and \c rhs.
  uqGslMatrixClass& operator+=(const uqGslMatrixClass& rhs);
  
  //! Stores in \c this the coordinate-wise subtraction of \c this by \c rhs.
  uqGslMatrixClass& operator-=(const uqGslMatrixClass& rhs);
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
  unsigned int      numRowsLocal              () const;
        
  //! Returns the global row dimension of \c this matrix.
  unsigned int      numRowsGlobal             () const;
        
  //! Returns the column dimension of \c this matrix.
  unsigned int      numCols                   () const;
  
  //! Returns the maximum element value of the matrix.
  double            max                       () const;
  
  //! This function returns the number of singular values of \c this matrix (rank). 
  /*! The rank function provides an estimate of the number of linearly independent rows or columns of a full matrix. */
  unsigned int      rank                      (double absoluteZeroThreshold, double relativeZeroThreshold) const;
  
     //! This function calculated the transpose of \c this matrix  (square).
  uqGslMatrixClass  transpose                 () const;
	
  //! This function calculated the inverse of \c this matrix (square).
  uqGslMatrixClass  inverse                   () const;
  
    
  //! Calculates the determinant of \c this matrix.
  double            determinant               () const;
	
  //! Calculates the ln(determinant) of \c this matrix.
  double            lnDeterminant             () const;
  
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

  //! Computes Cholesky factorization of a real symmetric positive definite matrix \c this. 
  /*! In case \this fails to be symmetric and positive definite, an error will be returned. */
  int               chol                      ();
	
//! Checks for the dimension of \c this matrix, \c matU, \c VecS and \c matVt, and calls the protected routine \c internalSvd to compute the singular values of \c this. 
  int               svd                       (uqGslMatrixClass& matU, uqGslVectorClass& vecS, uqGslMatrixClass& matVt) const;
        
  //! This function calls private member uqGslMatrixClass::internalSvd() to set a M-by-N orthogonal matrix U of the singular value decomposition (svd) of a general rectangular M-by-N matrix A. 
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  const uqGslMatrixClass& svdMatU                   () const;
  
  //! This function calls private member  uqGslMatrixClass::internalSvd() to set a N-by-N orthogonal square matrix V of the singular value decomposition (svd) of a general rectangular M-by-N matrix A. 
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  const uqGslMatrixClass& svdMatV                   () const;
  
  //! This function solves the system A x = b using the singular value decomposition (U, S, V) of A which must have been computed previously with uqGslMatrixClass::svd (x=solVec, b=rhsVec). 
  int               svdSolve                  (const uqGslVectorClass& rhsVec, uqGslVectorClass& solVec) const;
        
  //! This function solves the system A x = b using the singular value decomposition (U, S, V) of A which must have been computed previously with uqGslMatrixClass::svd (x=solMat, b=rhsMat).
  int               svdSolve                  (const uqGslMatrixClass& rhsMat, uqGslMatrixClass& solMat) const;        
  


  //! This function multiplies \c this matrix by vector \c x and returns the resulting vector.
  uqGslVectorClass  multiply                  (const uqGslVectorClass& x) const;

  //! This function calculates the inverse of \c this matrix and multiplies it with vector \c b. 
  /*! It calls void uqGslMatrixClass::invertMultiply(const uqGslVectorClass& b, uqGslVectorClass& x) internally.*/
  uqGslVectorClass  invertMultiply            (const uqGslVectorClass& b) const;
	
  //! This function calculates the inverse of \c this matrix, multiplies it with vector \c b and stores the result in vector \c x.
  /*! It checks for a previous LU decomposition of \c this matrix and does not recompute it
   if m_MU != NULL .*/
  void              invertMultiply            (const uqGslVectorClass& b, uqGslVectorClass& x) const;
	
  //! This function calculates the inverse of \c this matrix and multiplies it with matrix \c B.
  /*! It calls void uqGslMatrixClass::invertMultiply(const uqGslMatrixClass& B, uqGslMatrixClass& X) const internally.*/
  uqGslMatrixClass  invertMultiply            (const uqGslMatrixClass& B) const;

  //! This function calculates the inverse of \c this matrix, multiplies it with matrix \c B and stores the result in matrix \c X.
  /*! It checks for a previous LU decomposition of \c this matrix and does not recompute it
   if m_MU != NULL .*/
  void              invertMultiply            (const uqGslMatrixClass& B, uqGslMatrixClass& X) const;
  
  //! This function calculates the inverse of \c this matrix and multiplies it with vector \c b. 
  /*! It calls void uqGslMatrixClass::InvertMultiplyForceLU(const uqGslVectorClass& b, uqGslVectorClass& x) const;(const uqGslVectorClass& b,
uqGslVectorClass& x) internally.*/
  uqGslVectorClass  invertMultiplyForceLU     (const uqGslVectorClass& b) const;
	
  //! This function calculates the inverse of \c this matrix, multiplies it with vector \c b and stores the result in vector \c x.
  /*! It recalculates the LU decomposition of \c this matrix.*/	
  void              invertMultiplyForceLU     (const uqGslVectorClass& b, uqGslVectorClass& x) const;
  
    //! This function gets the column_num-th column of \c this matrix and stores it into vector \c column.
  void              getColumn                 (const unsigned int column_num, uqGslVectorClass& column) const;
  
  //! This function gets the column_num-th column of \c this matrix.
  uqGslVectorClass  getColumn                 (const unsigned int column_num) const;
  
  //! This function copies vector \c column into the column_num-th column of \c this matrix. 
  void              setColumn                 (const unsigned int column_num, const uqGslVectorClass& column);
  
  //! This function gets the row_num-th column of \c this matrix and stores it into vector \c row.
  void              getRow                    (const unsigned int row_num, uqGslVectorClass& row) const;
  
  //! This function gets the row_num-th column of \c this matrix.
  uqGslVectorClass  getRow                    (const unsigned int row_num) const;
  
  //! This function copies vector \c row into the row_num-th column of \c this matrix. 
  void              setRow                    (const unsigned int row_num, const uqGslVectorClass& row);
  
  //! This function computes the eigenvalues of a real symmetric matrix.
  void              eigen                     (uqGslVectorClass& eigenValues, uqGslMatrixClass* eigenVectors) const;
	
  //! This function finds largest eigenvalue, namely \c eigenValue, of \c this matrix and its corresponding eigenvector, namely \c eigenVector.
  void              largestEigen              (double& eigenValue, uqGslVectorClass& eigenVector) const;
 
  //! This function finds smallest eigenvalue, namely \c eigenValue, of \c this matrix and its corresponding eigenvector, namely \c eigenVector.
  void              smallestEigen             (double& eigenValue, uqGslVectorClass& eigenVector) const;

 
  //@} 
  
  //! @name Get/Set methods
  //@{  
  
  //! Component-wise set all values to \c this with value.
  void              cwSet                     (double value);
        
  //! Set the components of \c which positions are greater than (rowId,colId) with the value of mat(rowId,colId).
  void              cwSet                     (unsigned int rowId, unsigned int colId, const uqGslMatrixClass& mat);
  
  void              cwExtract                 (unsigned int rowId, unsigned int colId, uqGslMatrixClass& mat) const;
  
   //! This function sets all the entries bellow the main diagonal of \c this matrix to zero.
  /*! If \c includeDiagonal = false, then only the entries bellow the main diagonal are set to zero; 
  if \c includeDiagonal = true, then the elements of the matrix diagonal are also set to zero.*/
  void              zeroLower                 (bool includeDiagonal = false);
	
  //! This function sets all the entries above the main diagonal of \c this matrix to zero.
  /*! If \c includeDiagonal = false, then only the entries above the main diagonal are set to zero; 
  if \c includeDiagonal = true, then the elements of the matrix diagonal are also set to zero.*/
  void              zeroUpper                 (bool includeDiagonal = false);
  
  //! This function sets to zero (filters) all entries of \c this matrix which are smaller than \c thresholdValue.
  /*! If \c thresholdValue < 0 then no values will be filtered.*/
  void              filterSmallValues         (double thresholdValue);
  
  //! This function sets to zero (filters) all entries of \c this matrix which are greater than \c thresholdValue.
  /*! If \c thresholdValue < 0 then no values will be filtered.*/	
  void              filterLargeValues         (double thresholdValue);
  
  //! This function stores the transpose of \c this matrix into \c this matrix.
  void              fillWithTranspose         (unsigned int            rowId,
                                               unsigned int            colId,
                                               const uqGslMatrixClass& mat,
                                               bool                    checkForExactNumRowsMatching,
                                               bool                    checkForExactNumColsMatching);
  
  //! This function fills \c this matrix diagonally with const block  matrices.
  void              fillWithBlocksDiagonally  (unsigned int                                 rowId,
                                               unsigned int                                 colId,
                                               const std::vector<const uqGslMatrixClass* >& matrices,
                                               bool                                         checkForExactNumRowsMatching,
                                               bool                                         checkForExactNumColsMatching);
  
  //! This function fills \c this matrix diagonally with block matrices.
  void              fillWithBlocksDiagonally  (unsigned int                                 rowId,
                                               unsigned int                                 colId,
                                               const std::vector<      uqGslMatrixClass* >& matrices,
                                               bool                                         checkForExactNumRowsMatching,
                                               bool                                         checkForExactNumColsMatching);
  
  //! This function fills \c this matrix horizontally with const block matrices.
  void              fillWithBlocksHorizontally(unsigned int                                 rowId,
                                               unsigned int                                 colId,
                                               const std::vector<const uqGslMatrixClass* >& matrices,
                                               bool                                         checkForExactNumRowsMatching,
                                               bool                                         checkForExactNumColsMatching);
  
  //! This function fills \c this matrix horizontally with const block 
  void              fillWithBlocksHorizontally(unsigned int                                 rowId,
                                               unsigned int                                 colId,
                                               const std::vector<      uqGslMatrixClass* >& matrices,
                                               bool                                         checkForExactNumRowsMatching,
                                               bool                                         checkForExactNumColsMatching);
  
  //! This function fills \c this matrix vertically with const block matrices.
  void              fillWithBlocksVertically  (unsigned int                                 rowId,
                                               unsigned int                                 colId,
                                               const std::vector<const uqGslMatrixClass* >& matrices,
                                               bool                                         checkForExactNumRowsMatching,
                                               bool                                         checkForExactNumColsMatching);
  
  //! This function fills \c this matrix vertically with block matrices.
  void              fillWithBlocksVertically  (unsigned int                                 rowId,
                                               unsigned int                                 colId,
                                               const std::vector<      uqGslMatrixClass* >& matrices,
                                               bool                                         checkForExactNumRowsMatching,
                                               bool                                         checkForExactNumColsMatching);
	
  //! This function calculates the tensor product of matrices \c mat1 and \c mat2 and stores it in \c this matrix.
  void              fillWithTensorProduct     (unsigned int            rowId,
                                               unsigned int            colId,
                                               const uqGslMatrixClass& mat1,
                                               const uqGslMatrixClass& mat2,
                                               bool                    checkForExactNumRowsMatching,
                                               bool                    checkForExactNumColsMatching);

  //! This function calculates the tensor product of matrix \c mat1 and  vector \c vec2 and stores it in \c this matrix.
  void              fillWithTensorProduct     (unsigned int            rowId,
                                               unsigned int            colId,
                                               const uqGslMatrixClass& mat1,
                                               const uqGslVectorClass& vec2,
                                               bool                    checkForExactNumRowsMatching,
                                               bool                    checkForExactNumColsMatching);
  //@}	

 
  //! @name Miscellaneous methods
  //@{
  //! Returns \c this matrix.  
  gsl_matrix*       data                      ();
     
  void              mpiSum                    (const uqMpiCommClass& comm, uqGslMatrixClass& M_global) const;
  
  void              matlabLinearInterpExtrap  (const uqGslVectorClass& x1Vec, const uqGslMatrixClass& y1Mat, const uqGslVectorClass& x2Vec);
  //@}
	
        
	
  //! @name I/O methods
  //@{
    
  //! Print method. Defines the behavior of the ostream << operator inherited from the Object class.   
  void              print                     (std::ostream& os) const;
  
  //! Write contents of subenvironment in file \c fileName.
  void              subWriteContents          (const std::string&            varNamePrefix,
                                               const std::string&            fileName,
                                               const std::string&            fileType,
                                               const std::set<unsigned int>& allowedSubEnvIds) const;
					  
  //! Read contents of subenvironment from file \c fileName.
  void              subReadContents           (const std::string&            fileName,
                                               const std::string&            fileType,
                                               const std::set<unsigned int>& allowedSubEnvIds);
 //@}

private:
  //! In this function \c this matrix receives a copy of matrix \c src.
  void              copy                      (const uqGslMatrixClass& src);
  
  //! In this function resets the LU decomposition of \c this matrix, as well as deletes the private member pointers, if existing.
  void              resetLU                   ();
	
  //! This function multiplies \c this matrix by vector \c x and stores the resulting vector in \c y.
  void              multiply                  (const uqGslVectorClass& x, uqGslVectorClass& y) const;
        
  //! This function factorizes the M-by-N matrix A into the singular value decomposition A = U S V^T for M >= N. On output the matrix A is replaced by U.
  int               internalSvd               () const;

  //! GSL matrix, also referred to as \c this matrix.
          gsl_matrix*       m_mat;
	  
  //! GSL matrix for the LU decomposition of m_mat.	  
  mutable gsl_matrix*       m_LU;
  
  //! Inverse matrix of \c this.
  mutable uqGslMatrixClass* m_inverse;
  
  //! Mapping for matrices involved in the singular value decomposition (svd) routine.
  mutable uqMapClass*       m_svdColMap;
  
  //! m_svdUmat stores the M-by-N orthogonal matrix U after the singular value decomposition of a matrix.
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an 
   * M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S 
   * and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  mutable uqGslMatrixClass* m_svdUmat;
  
  //! m_svdSvec stores the diagonal of the N-by-N diagonal matrix of singular values S after the singular value decomposition of a matrix.
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an 
   * M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S 
   * and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  
  mutable uqGslVectorClass* m_svdSvec;
  
  //! m_svdVmat stores the N-by-N orthogonal square matrix V after the singular value decomposition of a matrix.
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an 
   * M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S 
   * and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  mutable uqGslMatrixClass* m_svdVmat;
  
  //! m_svdVmatT stores the transpose of N-by-N orthogonal square matrix V, namely V^T, after the singular value decomposition of a matrix.
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an 
   * M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S 
   * and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  mutable uqGslMatrixClass* m_svdVTmat;
  
  //! The determinant of \c this matrix.
  mutable double            m_determinant;
  
  //! The natural logarithm of the determinant of \c this matrix.
  mutable double            m_lnDeterminant;
  
  //! The permutation matrix of a LU decomposition.
  /*! In the  LU decomposition PA = LU, the j-th column of the matrix P is given by 
   * the k-th column of the identity matrix, where k = p_j the j-th element of the permutation vector. 
   * A is the square matrix of interest and P is the permutation matrix.     */
  mutable gsl_permutation*  m_permutation;
  
  //! m_signum stores the sign of the permutation of the LU decomposition PA = LU.
  /*! In the  LU decomposition PA = LU, where A is the square matrix of interest 
   * and P is the permutation matrix, m_signum has the value (-1)^n, 
   * where n is the number of interchanges in the permutation.*/
  mutable int               m_signum;
  
  //! Indicates whether or not \c this matrix is singular.
  mutable bool              m_isSingular;
};

uqGslMatrixClass operator*       (double a,                    const uqGslMatrixClass& mat);
uqGslVectorClass operator*       (const uqGslMatrixClass& mat, const uqGslVectorClass& vec);
uqGslMatrixClass operator*       (const uqGslMatrixClass& m1,  const uqGslMatrixClass& m2 );
uqGslMatrixClass operator+       (const uqGslMatrixClass& m1,  const uqGslMatrixClass& m2 );
uqGslMatrixClass operator-       (const uqGslMatrixClass& m1,  const uqGslMatrixClass& m2 );
uqGslMatrixClass matrixProduct   (const uqGslVectorClass& v1,  const uqGslVectorClass& v2 );

//! Row \c i of the returned matrix is equal to row \c i of \c mat multiplied by element \c i of \c v.
uqGslMatrixClass leftDiagScaling (const uqGslVectorClass& vec, const uqGslMatrixClass& mat);

//! Column \c j of the returned matrix is equal to column \c j of \c mat multiplied by element \c j of \c v.
uqGslMatrixClass rightDiagScaling(const uqGslMatrixClass& mat, const uqGslVectorClass& vec);
std::ostream&    operator<<      (std::ostream& os,            const uqGslMatrixClass& obj);

#endif // __UQ_GSL_MATRIX_H__
