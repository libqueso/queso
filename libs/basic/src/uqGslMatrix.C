/* libs/basic/src/uqGslMatrix.C
 * 
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <uqGslMatrix.h>
#include <uqGslVector.h>
#include <uqDefines.h>
#include <gsl/gsl_linalg.h>

uqGslMatrixClass::uqGslMatrixClass()
  :
  uqMatrixClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqGslMatrixClass::constructor(), default",
                      "should not be used by user");
}

uqGslMatrixClass::uqGslMatrixClass(
  const uqBaseEnvironmentClass& env,
  const Epetra_Map&         map,
  unsigned int              numCols)
  :
  uqMatrixClass(env,map),
  m_mat        (gsl_matrix_calloc(map.NumGlobalElements(),numCols)),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.rank(),
                      "uqGslMatrixClass::constructor()",
                      "null matrix generated");
}
 
uqGslMatrixClass::uqGslMatrixClass(
  const uqBaseEnvironmentClass& env,
  const Epetra_Map&         map,
  double                    diagValue)
  :
  uqMatrixClass(env,map),
  m_mat        (gsl_matrix_calloc(map.NumGlobalElements(),map.NumGlobalElements())),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.rank(),
                      "uqGslMatrixClass::constructor(), eye",
                      "null matrix generated");

  for (unsigned int i = 0; i < m_mat->size1; ++i) {
    (*this)(i,i) = diagValue;
  }
}

uqGslMatrixClass::uqGslMatrixClass(
  const uqGslVectorClass& v,
  double                  diagValue)
  :
  uqMatrixClass(v.env(),v.map()),
  m_mat        (gsl_matrix_calloc(v.size(),v.size())),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.rank(),
                      "uqGslMatrixClass::constructor(), eye",
                      "null matrix generated");

  for (unsigned int i = 0; i < m_mat->size1; ++i) {
    (*this)(i,i) = diagValue;
  }
}

uqGslMatrixClass::uqGslMatrixClass(const uqGslVectorClass& v)
  :
  uqMatrixClass(v.env(),v.map()),
  m_mat        (gsl_matrix_calloc(v.size(),v.size())),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.rank(),
                      "uqGslMatrixClass::constructor(), from vector",
                      "null matrix generated");

  unsigned int dim = v.size();
  for (unsigned int i = 0; i < dim; ++i) {
    (*this)(i,i) = v[i];
  }
}

uqGslMatrixClass::uqGslMatrixClass(const uqGslMatrixClass& B)
  :
  uqMatrixClass(B.env(),B.map()),
  m_mat        (gsl_matrix_calloc(B.numRows(),B.numCols())),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.rank(),
                      "uqGslMatrixClass::constructor(), copy",
                      "null vector generated");
  this->uqMatrixClass::copy(B);
  this->copy(B);
}

uqGslMatrixClass::~uqGslMatrixClass()
{
  if (m_LU)          gsl_matrix_free     (m_LU);
  if (m_permutation) gsl_permutation_free(m_permutation);
  if (m_mat)         gsl_matrix_free     (m_mat);
}

uqGslMatrixClass&
uqGslMatrixClass::operator=(const uqGslMatrixClass& obj)
{
  this->copy(obj);
  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator*=(double a)
{
  int iRC;
  iRC = gsl_matrix_scale(m_mat,a);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqGslMatrixClass::operator*=()",
                    "scaling failed");
  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator/=(double a)
{
  *this *= (1./a);

  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator+=(const uqGslMatrixClass& rhs)
{
  int iRC;
  iRC = gsl_matrix_add(m_mat,rhs.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqGslMatrixClass::operator+=()",
                    "failed");

  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator-=(const uqGslMatrixClass& rhs)
{
  int iRC;
  iRC = gsl_matrix_sub(m_mat,rhs.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqGslMatrixClass::operator-=()",
                    "failed");

  return *this;
}

double&
uqGslMatrixClass::operator()(unsigned int i, unsigned int j)
{
  return *gsl_matrix_ptr(m_mat,i,j);
}

const double&
uqGslMatrixClass::operator()(unsigned int i, unsigned int j) const
{
  return *gsl_matrix_const_ptr(m_mat,i,j);
}

void
uqGslMatrixClass::copy(const uqGslMatrixClass& src)
{
  int iRC;
  iRC = gsl_matrix_memcpy(this->m_mat, src.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqGslMatrixClass::copy()",
                    "failed");

  return;
}

unsigned int
uqGslMatrixClass::numRows() const
{
  return m_mat->size1;
}

unsigned int
uqGslMatrixClass::numCols() const
{
  return m_mat->size2;
}

int
uqGslMatrixClass::chol()
{
  int iRC;
  //std::cout << "Calling gsl_linalg_cholesky_decomp()..." << std::endl;
  gsl_error_handler_t* oldHandler;
  oldHandler = gsl_set_error_handler_off();
  iRC = gsl_linalg_cholesky_decomp(m_mat);
  gsl_set_error_handler(oldHandler);
  //std::cout << "Returned from gsl_linalg_cholesky_decomp() with iRC = " << iRC << std::endl;
  UQ_RC_MACRO(iRC, // Yes, *not* a fatal check on RC
              m_env.rank(),
              "uqGslMatrixClass::chol()",
              "matrix is not positive definite",
              UQ_MATRIX_IS_NOT_POS_DEFINITE_RC);

  return iRC;
}

void
uqGslMatrixClass::zeroLower(bool includeDiagonal)
{
  if (this->numRows() != this->numCols()) return;

  unsigned int dim = this->numRows();
  if (includeDiagonal) {
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j <= i; j++) {
        (*this)(i,j) = 0.;
      }
    }
  }
  else {
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = 0; j < i; j++) {
        (*this)(i,j) = 0.;
      }
    }
  }

  return;
}

void
uqGslMatrixClass::zeroUpper(bool includeDiagonal)
{
  if (this->numRows() != this->numCols()) return;

  unsigned int dim = this->numRows();
  if (includeDiagonal) {
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = i; j < dim; j++) {
        (*this)(i,j) = 0.;
      }
    }
  }
  else {
    for (unsigned int i = 0; i < dim; i++) {
      for (unsigned int j = (i+1); j < dim; j++) {
        (*this)(i,j) = 0.;
      }
    }
  }

  return;
}

uqGslMatrixClass
uqGslMatrixClass::transpose () const
{
  unsigned int nRows = this->numRows();
  unsigned int nCols = this->numCols();

  uqGslMatrixClass mat(m_env,m_map,nCols);
  for (unsigned int row = 0; row < nRows; ++row) {
    for (unsigned int col = 0; col < nCols; ++col) {
      mat(row,col) = (*this)(col,row);
    }
  }

  return mat;
}

uqGslVectorClass
uqGslMatrixClass::multiply(
  const uqGslVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != x.size()),
                      m_env.rank(),
                      "uqGslMatrixClass::multiply(), return vector",
                      "matrix and vector have incompatible sizes");

  uqGslVectorClass y(m_env,m_map);
  this->multiply(x,y);

  return y;
}

void
uqGslMatrixClass::multiply(
  const uqGslVectorClass& x,
        uqGslVectorClass& y) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != x.size()),
                      m_env.rank(),
                      "uqGslMatrixClass::multiply(), vector return void",
                      "matrix and x have incompatible sizes");

  UQ_FATAL_TEST_MACRO((y.size() != x.size()),
                      m_env.rank(),
                      "uqGslMatrixClass::multiply(), vector return void",
                      "y and x have incompatible sizes");

  unsigned int sizeX = this->numCols();
  unsigned int sizeY = this->numRows();
  for (unsigned int i = 0; i < sizeY; ++i) {
    double value = 0.;
    for (unsigned int j = 0; j < sizeX; ++j) {
      value += (*this)(i,j)*x[j];
    }
    y[i] = value;
  }

  return;
}

uqGslVectorClass
uqGslMatrixClass::invertMultiply(
  const uqGslVectorClass& b) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.size()),
                      m_env.rank(),
                      "uqGslMatrixClass::invertMultiply(), return vector",
                      "matrix and rhs have incompatible sizes");

  uqGslVectorClass x(m_env,m_map);
  this->invertMultiply(b,x);

  return x;
}

void
uqGslMatrixClass::invertMultiply(
  const uqGslVectorClass& b,
        uqGslVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.size()),
                      m_env.rank(),
                      "uqGslMatrixClass::multiply(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.size() != b.size()),
                      m_env.rank(),
                      "uqGslMatrixClass::multiply(), return void",
                      "solution and rhs have incompatible sizes");

  int iRC;
  if (m_LU == NULL) {
    m_LU = gsl_matrix_calloc(numRows(),numCols());
    UQ_FATAL_TEST_MACRO((m_LU == NULL),
                        m_env.rank(),
                        "uqGslMatrixClass::invertMultiply()",
                        "gsl_matrix_calloc() failed");

    iRC = gsl_matrix_memcpy(m_LU, m_mat);
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.rank(),
                      "uqGslMatrixClass::invertMultiply()",
                      "gsl_matrix_memcpy() failed");

    m_permutation = gsl_permutation_calloc(numCols());
    UQ_FATAL_TEST_MACRO((m_permutation == NULL),
                        m_env.rank(),
                        "uqGslMatrixClass::invertMultiply()",
                        "gsl_permutation_calloc() failed");

    int signum;
    iRC = gsl_linalg_LU_decomp(m_LU,m_permutation,&signum); 
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.rank(),
                      "uqGslMatrixClass::invertMultiply()",
                      "gsl_linalg_LU_decomp() failed");
  }

  iRC = gsl_linalg_LU_solve(m_LU,m_permutation,b.data(),x.data()); 
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.rank(),
                    "uqGslMatrixClass::invertMultiply()",
                    "gsl_linalg_LU_solve() failed");

  return;
}

void
uqGslMatrixClass::print(std::ostream& os) const
{
  unsigned int nRows = this->numRows();
  unsigned int nCols = this->numCols();

  if (m_printHorizontally) {
    for (unsigned int i = 0; i < nRows; ++i) {
      for (unsigned int j = 0; j < nCols; ++j) {
        os << (*this)(i,j)
           << " ";
      }
      os << "# ";
    }
    //os << std::endl;
  }
  else {
    for (unsigned int i = 0; i < nRows; ++i) {
      for (unsigned int j = 0; j < nCols; ++j) {
        os << (*this)(i,j)
           << " ";
      }
      os << std::endl;
    }
  }

  return;
}

std::ostream&
operator<<(std::ostream& os, const uqGslMatrixClass& obj)
{
  obj.print(os);

  return os;
}

uqGslMatrixClass operator*(double a, const uqGslMatrixClass& mat)
{
  uqGslMatrixClass answer(mat);
  answer *= a;
  return answer;
}

uqGslVectorClass operator*(const uqGslMatrixClass& mat, const uqGslVectorClass& vec)
{
  return mat.multiply(vec);
}

uqGslMatrixClass operator*(const uqGslMatrixClass& m1, const uqGslMatrixClass& m2)
{
  unsigned int m1Rows = m1.numRows();
  unsigned int m1Cols = m1.numCols();
  unsigned int m2Rows = m2.numRows();
  unsigned int m2Cols = m2.numCols();

  UQ_FATAL_TEST_MACRO((m1Cols != m2Rows),
                      m1.env().rank(),
                      "uqGslMatrixClass operator*(matrix,matrix)",
                      "different sizes m1Cols and m2Rows");

  uqGslMatrixClass mat(m1.env(),m1.map(),m2Cols);

  unsigned int commonSize = m1Cols;
  for (unsigned int row1 = 0; row1 < m1Rows; ++row1) {
    for (unsigned int col2 = 0; col2 < m2Cols; ++col2) {
      double result = 0.;
      for (unsigned int k = 0; k < commonSize; ++k) {
        result += m1(row1,k)*m2(k,col2);
      }
      mat(row1,col2) = result;
    }
  }

  return mat;
}

uqGslMatrixClass operator+(const uqGslMatrixClass& m1, const uqGslMatrixClass& m2)
{
  uqGslMatrixClass answer(m1);
  answer += m2;
  return answer;
}

uqGslMatrixClass matrixProduct(const uqGslVectorClass& v1, const uqGslVectorClass& v2)
{
  unsigned int numRows = v1.size();
  unsigned int numCols = v2.size();
  uqGslMatrixClass answer(v1.env(),v1.map(),numCols);

  for (unsigned int i = 0; i < numRows; ++i) {
    double value1 = v1[i];
    for (unsigned int j = 0; j < numCols; ++j) {
      answer(i,j) = value1*v2[j];
    }
  }

  return answer;
}

uqGslMatrixClass diagScaling(const uqGslVectorClass& vec, const uqGslMatrixClass& mat)
{
  unsigned int vSize = vec.size();
  unsigned int mRows = mat.numRows();
  unsigned int mCols = mat.numCols();

  UQ_FATAL_TEST_MACRO((vSize != mRows),
                      mat.env().rank(),
                      "uqGslMatrixClass diagScaling(vector,matrix)",
                      "size of vector is different from the number of rows in matrix");

  uqGslMatrixClass answer(mat);
  for (unsigned int i = 0; i < mRows; ++i) {
    double vecValue = vec[i];
    for (unsigned int j = 0; j < mCols; ++j) {
      answer(i,j) *= vecValue;
    }
  }

  return answer;
}
