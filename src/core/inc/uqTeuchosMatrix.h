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

class uqTeuchosMatrixClass : public uqMatrixClass
{
public:
  uqTeuchosMatrixClass();
  uqTeuchosMatrixClass(const uqBaseEnvironmentClass& env,
                       const uqMapClass&             map,
                       unsigned int                  numCols);
  uqTeuchosMatrixClass(const uqBaseEnvironmentClass& env,
                       const uqMapClass&             map,
     	               double                        diagValue); 
  uqTeuchosMatrixClass(const uqTeuchosVectorClass&   v,
                   	   double                        diagValue); 
  uqTeuchosMatrixClass(const uqTeuchosVectorClass&   v);         
  uqTeuchosMatrixClass(const uqTeuchosMatrixClass&   B);
 ~uqTeuchosMatrixClass();


  uqTeuchosMatrixClass& operator= (const uqTeuchosMatrixClass& rhs);
  uqTeuchosMatrixClass& operator*=(double a);
  uqTeuchosMatrixClass& operator/=(double a);
  uqTeuchosMatrixClass& operator+=(const uqTeuchosMatrixClass& rhs);
  uqTeuchosMatrixClass& operator-=(const uqTeuchosMatrixClass& rhs);
            double& operator()(unsigned int i, unsigned int j);
      const double& operator()(unsigned int i, unsigned int j) const;


  unsigned int      numRowsLocal           () const;
  unsigned int      numRowsGlobal          () const;
  unsigned int      numCols                () const;
  double*           values                 () ;//added by Kemelli on 12/04/12
  int 	            stride                 () ;//added by Kemelli on 12/04/12
  int               chol                   () ; 
  void              zeroLower              (bool includeDiagonal = false);
  void              zeroUpper              (bool includeDiagonal = false);
  void              print                  (std::ostream& os) const;
  void              mpiSum                 (const uqMpiCommClass& comm, uqTeuchosMatrixClass& M_global) const;
  int               svd                    (uqTeuchosMatrixClass& matU, uqTeuchosVectorClass& vecS, uqTeuchosMatrixClass& matVt) const;
  const uqTeuchosMatrixClass& svdMatU      () const;
  const uqTeuchosMatrixClass& svdMatV      () const;
  int               svdSolve               (const uqTeuchosVectorClass& rhsVec, uqTeuchosVectorClass& solVec) const;
  int               svdSolve               (const uqTeuchosMatrixClass& rhsMat, uqTeuchosMatrixClass& solMat) const;

  uqTeuchosMatrixClass  transpose                 () const;
  uqTeuchosMatrixClass  inverse                   () const;
  uqTeuchosVectorClass  multiply                  (const uqTeuchosVectorClass& x) const;
  uqTeuchosMatrixClass  invertMultiply            (const uqTeuchosMatrixClass& B) const;
  void                  invertMultiply            (const uqTeuchosMatrixClass& B, uqTeuchosMatrixClass& X) const;
  uqTeuchosVectorClass  invertMultiply            (const uqTeuchosVectorClass& b) const;
  void                  invertMultiply            (const uqTeuchosVectorClass& b, uqTeuchosVectorClass& x) const;
  uqTeuchosVectorClass  invertMultiplyForceLU     (const uqTeuchosVectorClass& b) const;
  void                  invertMultiplyForceLU     (const uqTeuchosVectorClass& b, uqTeuchosVectorClass& x) const;
  double                determinant               () const;
  double                lnDeterminant             () const;
  void                  getColumn                 (const unsigned int column_num, uqTeuchosVectorClass& column) const;
  uqTeuchosVectorClass  getColumn                 (const unsigned int column_num) const;
  void                  setColumn                 (const unsigned int column_num, const uqTeuchosVectorClass& column);
  void                  getRow                    (const unsigned int row_num, uqTeuchosVectorClass& row) const;
  uqTeuchosVectorClass  getRow                    (const unsigned int row_num) const;
  void                  setRow                    (const unsigned int row_num, const uqTeuchosVectorClass& row);
  void                  subReadContents           (const std::string&            fileName,
                                                   const std::string&            fileType,
                                                   const std::set<unsigned int>& allowedSubEnvIds);
  void                  subWriteContents          (const std::string&            varNamePrefix,
                                                   const std::string&            fileName,
                                                   const std::string&            fileType,
                                                   const std::set<unsigned int>& allowedSubEnvIds) const;
						     
  void              eigen                     (uqTeuchosVectorClass& eigenValues, uqTeuchosMatrixClass* eigenVectors) const;			
  void              largestEigen              (double& eigenValue, uqTeuchosVectorClass& eigenVector) const;
  void              smallestEigen             (double& eigenValue, uqTeuchosVectorClass& eigenVector) const;
  
  double            normFrob                  () const;
  double            normMax                   () const;
  double            max                       () const;
  unsigned int      rank                      (double absoluteZeroThreshold, double relativeZeroThreshold) const;
  void              cwSet                     (double value);
  void              cwSet                     (unsigned int rowId, unsigned int colId, const uqTeuchosMatrixClass& mat);

  void              filterSmallValues         (double thresholdValue);
  void              filterLargeValues         (double thresholdValue);
  void              fillWithTranspose         (const uqTeuchosMatrixClass& mat); 
  void              fillWithBlocksDiagonally  (const std::vector<const uqTeuchosMatrixClass* >& matrices);
  void              fillWithBlocksDiagonally  (const std::vector<      uqTeuchosMatrixClass* >& matrices);
  void              fillWithBlocksHorizontally(const std::vector<const uqTeuchosMatrixClass* >& matrices);
  void              fillWithBlocksHorizontally(const std::vector<      uqTeuchosMatrixClass* >& matrices);
  void              fillWithBlocksVertically  (const std::vector<const uqTeuchosMatrixClass* >& matrices);
  void              fillWithBlocksVertically  (const std::vector<      uqTeuchosMatrixClass* >& matrices);
  void              fillWithTensorProduct     (const uqTeuchosMatrixClass& mat1, const uqTeuchosMatrixClass& mat2);
  void              fillWithTensorProduct     (const uqTeuchosMatrixClass& mat1, const uqTeuchosVectorClass& vec2); 
  
#if 0
      void              matlabLinearInterpExtrap  (const uqTeuchosVectorClass& x1Vec, const uqTeuchosMatrixClass& y1Mat, const uqTeuchosVectorClass& x2Vec);

        Teuchos::SerialDenseMatrix<int,double>       data                      ();
#endif
private:
        void              copy                      (const uqTeuchosMatrixClass& src);
        void              resetLU                   ();
        void              multiply                  (const uqTeuchosVectorClass& x, uqTeuchosVectorClass& y) const;
        int               internalSvd               () const;

  Teuchos::SerialDenseMatrix<int,double> m_mat;
  mutable Teuchos::SerialDenseMatrix<int,double> m_LU; 
  //mutable uqTeuchosMatrixClass*	m_LU; 
  
  mutable uqTeuchosMatrixClass* m_inverse;
  mutable uqMapClass*       	m_svdColMap;
  mutable uqTeuchosMatrixClass* m_svdUmat;
  mutable uqTeuchosVectorClass* m_svdSvec;
  mutable uqTeuchosMatrixClass* m_svdVmat;
  mutable uqTeuchosMatrixClass* m_svdVTmat;
  mutable double                m_determinant;
  mutable double                m_lnDeterminant;
  mutable uqTeuchosMatrixClass*  m_permutation;
  mutable int*              v_pivoting; //Kemelli, added on 12/07/12
  mutable int               m_signum;
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

