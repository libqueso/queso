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
    \brief This class creates and provides basic support for matrices of templated 
    type as a specialization of uqMatrixClass using GSL matrices, which are defined 
    by an encapsulated gsl_matrix structure.
*/

class uqGslMatrixClass : public uqMatrixClass
{
public:
 //! @name Constructor/Destructor methods.
  //@{ 

  //! Default Constructor
  /*! Creates an empty matrix vector of no dimension. The Shaping methods should be used to size this matrix.
    Values of this matrix should be set using the [], (), or = operators.*/
  uqGslMatrixClass();
  
   //! Shaped Constructor
  
  //!Creates a shaped matrix with \c numCols cols.  All values are initialized to 0 when \c zeroOut is true.
  /*!   \param numCols - Number of columns in matrix. */  
  uqGslMatrixClass(const uqBaseEnvironmentClass& env,
                   const uqMapClass&             map,
                   unsigned int                  numCols);
  uqGslMatrixClass(const uqBaseEnvironmentClass& env,
                   const uqMapClass&             map,
                   double                        diagValue); // MATLAB eye
  uqGslMatrixClass(const uqGslVectorClass&       v,
                   double                        diagValue); // MATLAB eye
  uqGslMatrixClass(const uqGslVectorClass&       v);         // MATLAB diag
  uqGslMatrixClass(const uqGslMatrixClass&       B);
  
  //! Destructor
  ~uqGslMatrixClass();
  //@}


  //! @name Set methods.
  //@{ 
  //! 	Copies values from matrix rhs to \c this. 
  uqGslMatrixClass& operator= (const uqGslMatrixClass& rhs);
  
  //! Stores in \c this the cordinate-wise multiplication of \c this and a.
  uqGslMatrixClass& operator*=(double a);
  
  //! Stores in \c this the cordinate-wise division of \c this by a.
  uqGslMatrixClass& operator/=(double a);
  
  //! Stores in \c this the cordinate-wise addition of \c this and rhs.
  uqGslMatrixClass& operator+=(const uqGslMatrixClass& rhs);
  
  //! Stores in \c this the cordinate-wise subtraction of \c this by rhs.
  uqGslMatrixClass& operator-=(const uqGslMatrixClass& rhs);
  //@}
  
  
  //! @name Accessor methods.
  //@{
  //! Element access method (non-const).
            double& operator()(unsigned int i, unsigned int j);
            
  //! Element access method (const).  
        const double& operator()(unsigned int i, unsigned int j) const;
  //@}

  //! @name Attribute methods.
  //@{ 
  //! Returns the local row dimension of this matrix.
        unsigned int      numRowsLocal              () const;
        
  //! Returns the global row dimension of this matrix.
        unsigned int      numRowsGlobal             () const;
        
  //! Returns the column dimension of this matrix.
        unsigned int      numCols                   () const;
  //@}
  
  //! @name Norm methods.
  //@{ 
  //! Returns the Frobenius norm of the matrix.
        double            normFrob                  () const;
  //! Returns the Frobenius norm of the matrix.
        double            normMax                   () const;
  //@}
       

  
  //! @name Set methods.
  //@{ 
  //! Component-wise set all values to \c this with value.
        void              cwSet                     (double value);
        
  //! Set the componts of \c which positions are greater than (rowId,colId) with the value of mat(rowId,colId).
  // TODO: improve this description
        void              cwSet                     (unsigned int rowId, unsigned int colId, const uqGslMatrixClass& mat);
  //@}
  
  
  //! @name Mathematical methods.
  //@{ 
  //! Returns the maximum element value of the matrix.
        double            max                       () const;
  
  //! Cholesky decomposition of \thic as long as \this is a symmetric, positive definite square matrix
        int               chol                      ();
        int               svd                       (uqGslMatrixClass& matU, uqGslVectorClass& vecS, uqGslMatrixClass& matVt) const;
        
  //! This function solves the system A x = b using the singular value decomposition (U, S, V) of A which must have been computed previously with uqGslMatrixClass::svd (x=solVec, b=rhsVec). 
        int               svdSolve                  (const uqGslVectorClass& rhsVec, uqGslVectorClass& solVec) const;
        
  //! This function solves the system A x = b using the singular value decomposition (U, S, V) of A which must have been computed previously with uqGslMatrixClass::svd (x=solMat, b=rhsMat).
        int               svdSolve                  (const uqGslMatrixClass& rhsMat, uqGslMatrixClass& solMat) const;
        
        
  //! This function calls private member uqGslMatrixClass::internalSvd() to set a M-by-N orthogonal matrix U of the singular value decomposition (svd) of a general rectangular M-by-N matrix A. 
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  const uqGslMatrixClass& svdMatU                   () const;
  
  //! This function calls private member  uqGslMatrixClass::internalSvd() to set a N-by-N orthogonal square matrix V of the singular value decomposition (svd) of a general rectangular M-by-N matrix A. 
  /*! A general rectangular M-by-N matrix A has a singular value decomposition (svd) into the product of an M-by-N orthogonal matrix U, an N-by-N diagonal matrix of singular values S and the transpose of an N-by-N orthogonal square matrix V,  A = U S V^T.  */
  const uqGslMatrixClass& svdMatV                   () const;
  
        void              zeroLower                 (bool includeDiagonal = false);
        void              zeroUpper                 (bool includeDiagonal = false);
        
  //@}      
        void              filterSmallValues         (double thresholdValue);
        void              filterLargeValues         (double thresholdValue);
        uqGslMatrixClass  transpose                 () const;
        uqGslMatrixClass  inverse                   () const;
        void              fillWithBlocksDiagonally  (const std::vector<const uqGslMatrixClass* >& matrices);
        void              fillWithBlocksDiagonally  (const std::vector<      uqGslMatrixClass* >& matrices);
        void              fillWithBlocksHorizontally(const std::vector<const uqGslMatrixClass* >& matrices);
        void              fillWithBlocksHorizontally(const std::vector<      uqGslMatrixClass* >& matrices);
        void              fillWithBlocksVertically  (const std::vector<const uqGslMatrixClass* >& matrices);
        void              fillWithBlocksVertically  (const std::vector<      uqGslMatrixClass* >& matrices);
        void              fillWithTensorProduct     (const uqGslMatrixClass& mat1, const uqGslMatrixClass& mat2);
        void              fillWithTensorProduct     (const uqGslMatrixClass& mat1, const uqGslVectorClass& vec2);
        void              fillWithTranspose         (const uqGslMatrixClass& mat);
        double            determinant               () const;
        double            lnDeterminant             () const;
        unsigned int      rank                      (double absoluteZeroThreshold, double relativeZeroThreshold) const;
        uqGslVectorClass  multiply                  (const uqGslVectorClass& x) const;
        uqGslVectorClass  invertMultiply            (const uqGslVectorClass& b) const;
        void              invertMultiply            (const uqGslVectorClass& b, uqGslVectorClass& x) const;
        uqGslMatrixClass  invertMultiply            (const uqGslMatrixClass& B) const;
        void              invertMultiply            (const uqGslMatrixClass& B, uqGslMatrixClass& X) const;

        void              eigen                     (uqGslVectorClass& eigenValues, uqGslMatrixClass* eigenVectors) const;
        void              largestEigen              (double& eigenValue, uqGslVectorClass& eigenVector) const;
        void              smallestEigen             (double& eigenValue, uqGslVectorClass& eigenVector) const;

        uqGslVectorClass  invertMultiplyForceLU     (const uqGslVectorClass& b) const;
        void              invertMultiplyForceLU     (const uqGslVectorClass& b, uqGslVectorClass& x) const;

        void              print                     (std::ostream& os) const;
        void              subWriteContents          (const std::string&            varNamePrefix,
                                                     const std::string&            fileName,
                                                     const std::string&            fileType,
                                                     const std::set<unsigned int>& allowedSubEnvIds) const;
        void              subReadContents           (const std::string&            fileName,
                                                     const std::string&            fileType,
                                                     const std::set<unsigned int>& allowedSubEnvIds);

        void              getColumn                 (const unsigned int column_num, uqGslVectorClass& column) const;
        void              getRow                    (const unsigned int row_num, uqGslVectorClass& row) const;
        uqGslVectorClass  getColumn                 (const unsigned int column_num) const;
        uqGslVectorClass  getRow                    (const unsigned int row_num) const;
        void              setColumn                 (const unsigned int column_num, const uqGslVectorClass& column);
        void              setRow                    (const unsigned int row_num, const uqGslVectorClass& row);

        void              mpiSum                    (const uqMpiCommClass& comm, uqGslMatrixClass& M_global) const;
        void              matlabLinearInterpExtrap  (const uqGslVectorClass& x1Vec, const uqGslMatrixClass& y1Mat, const uqGslVectorClass& x2Vec);

        gsl_matrix*       data                      ();

private:
        void              copy                      (const uqGslMatrixClass& src);
        void              resetLU                   ();
        void              multiply                  (const uqGslVectorClass& x, uqGslVectorClass& y) const;
        
        //! This function factorizes the M-by-N matrix A into the singular value decomposition A = U S V^T for M >= N. On output the matrix A is replaced by U.
        int               internalSvd               () const;

          gsl_matrix*       m_mat;
  mutable gsl_matrix*       m_LU;
  mutable uqGslMatrixClass* m_inverse;
  mutable uqMapClass*       m_svdColMap;
  mutable uqGslMatrixClass* m_svdUmat;
  mutable uqGslVectorClass* m_svdSvec;
  mutable uqGslMatrixClass* m_svdVmat;
  mutable uqGslMatrixClass* m_svdVTmat;
  mutable double            m_determinant;
  mutable double            m_lnDeterminant;
  mutable gsl_permutation*  m_permutation;
  mutable int               m_signum;
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
