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

#include <uqMatrix.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <uqGslVector.h>

class uqGslMatrixClass : public uqMatrixClass
{
public:
  uqGslMatrixClass();
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
 ~uqGslMatrixClass();

  uqGslMatrixClass& operator= (const uqGslMatrixClass& rhs);
  uqGslMatrixClass& operator*=(double a);
  uqGslMatrixClass& operator/=(double a);
  uqGslMatrixClass& operator+=(const uqGslMatrixClass& rhs);
  uqGslMatrixClass& operator-=(const uqGslMatrixClass& rhs);
            double& operator()(unsigned int i, unsigned int j);
      const double& operator()(unsigned int i, unsigned int j) const;

        unsigned int      numRowsLocal              () const;
        unsigned int      numRowsGlobal             () const;
        unsigned int      numCols                   () const;
        double            normFrob                  () const;
        double            normMax                   () const;
        double            max                       () const;
        void              cwSet                     (double value);
        void              cwSet                     (unsigned int rowId, unsigned int colId, const uqGslMatrixClass& mat);
        int               chol                      ();
        int               svd                       (uqGslMatrixClass& matU, uqGslVectorClass& vecS, uqGslMatrixClass& matVt) const;
        int               svdSolve                  (const uqGslVectorClass& rhsVec, uqGslVectorClass& solVec) const;
        int               svdSolve                  (const uqGslMatrixClass& rhsMat, uqGslMatrixClass& solMat) const;
  const uqGslMatrixClass& svdMatU                   () const;
  const uqGslMatrixClass& svdMatV                   () const;
        void              zeroLower                 (bool includeDiagonal = false);
        void              zeroUpper                 (bool includeDiagonal = false);
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
        unsigned int      rank                      (double zeroThreshold) const;
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
};

uqGslMatrixClass operator*       (double a,                    const uqGslMatrixClass& mat);
uqGslVectorClass operator*       (const uqGslMatrixClass& mat, const uqGslVectorClass& vec);
uqGslMatrixClass operator*       (const uqGslMatrixClass& m1,  const uqGslMatrixClass& m2 );
uqGslMatrixClass operator+       (const uqGslMatrixClass& m1,  const uqGslMatrixClass& m2 );
uqGslMatrixClass operator-       (const uqGslMatrixClass& m1,  const uqGslMatrixClass& m2 );
uqGslMatrixClass matrixProduct   (const uqGslVectorClass& v1,  const uqGslVectorClass& v2 );
uqGslMatrixClass leftDiagScaling (const uqGslVectorClass& vec, const uqGslMatrixClass& mat);
uqGslMatrixClass rightDiagScaling(const uqGslMatrixClass& mat, const uqGslVectorClass& vec);
std::ostream&    operator<<      (std::ostream& os,            const uqGslMatrixClass& obj);

#endif // __UQ_GSL_MATRIX_H__
