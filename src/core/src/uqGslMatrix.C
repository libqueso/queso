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

#include <uqGslMatrix.h>
#include <uqGslVector.h>
#include <uqDefines.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <sys/time.h>
#include <cmath>

uqGslMatrixClass::uqGslMatrixClass()
  :
  uqMatrixClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqGslMatrixClass::constructor(), default",
                      "should not be used by user");
}

uqGslMatrixClass::uqGslMatrixClass( // can be a rectangular matrix
  const uqBaseEnvironmentClass& env,
  const uqMapClass&             map,
  unsigned int                  nCols)
  :
  uqMatrixClass  (env,map),
  m_mat          (gsl_matrix_calloc(map.NumGlobalElements(),nCols)),
  m_LU           (NULL),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  m_permutation  (NULL),
  m_signum       (0)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.worldRank(),
                      "uqGslMatrixClass::constructor()",
                      "null matrix generated");
}
 
uqGslMatrixClass::uqGslMatrixClass( // square matrix
  const uqBaseEnvironmentClass& env,
  const uqMapClass&             map,
  double                        diagValue)
  :
  uqMatrixClass  (env,map),
  m_mat          (gsl_matrix_calloc(map.NumGlobalElements(),map.NumGlobalElements())),
  m_LU           (NULL),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  m_permutation  (NULL),
  m_signum       (0)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.worldRank(),
                      "uqGslMatrixClass::constructor(), eye",
                      "null matrix generated");

  for (unsigned int i = 0; i < m_mat->size1; ++i) {
    (*this)(i,i) = diagValue;
  }
}

uqGslMatrixClass::uqGslMatrixClass( // square matrix
  const uqGslVectorClass& v,
  double                  diagValue)
  :
  uqMatrixClass  (v.env(),v.map()),
  m_mat          (gsl_matrix_calloc(v.sizeLocal(),v.sizeLocal())),
  m_LU           (NULL),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  m_permutation  (NULL),
  m_signum       (0)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.worldRank(),
                      "uqGslMatrixClass::constructor(), eye",
                      "null matrix generated");

  for (unsigned int i = 0; i < m_mat->size1; ++i) {
    (*this)(i,i) = diagValue;
  }
}

uqGslMatrixClass::uqGslMatrixClass(const uqGslVectorClass& v) // square matrix
  :
  uqMatrixClass  (v.env(),v.map()),
  m_mat          (gsl_matrix_calloc(v.sizeLocal(),v.sizeLocal())),
  m_LU           (NULL),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  m_permutation  (NULL),
  m_signum       (0)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.worldRank(),
                      "uqGslMatrixClass::constructor(), from vector",
                      "null matrix generated");

  unsigned int dim = v.sizeLocal();
  for (unsigned int i = 0; i < dim; ++i) {
    (*this)(i,i) = v[i];
  }
}

uqGslMatrixClass::uqGslMatrixClass(const uqGslMatrixClass& B) // can be a rectangular matrix
  :
  uqMatrixClass  (B.env(),B.map()),
  m_mat          (gsl_matrix_calloc(B.numRowsLocal(),B.numCols())),
  m_LU           (NULL),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  m_permutation  (NULL),
  m_signum       (0)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.worldRank(),
                      "uqGslMatrixClass::constructor(), copy",
                      "null vector generated");
  this->uqMatrixClass::copy(B);
  this->copy(B);
}

uqGslMatrixClass::~uqGslMatrixClass()
{
  this->resetLU();
  if (m_mat) gsl_matrix_free(m_mat);
}

uqGslMatrixClass&
uqGslMatrixClass::operator=(const uqGslMatrixClass& obj)
{
  this->resetLU();
  this->copy(obj);
  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator*=(double a)
{
  this->resetLU();
  int iRC;
  iRC = gsl_matrix_scale(m_mat,a);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
                    "uqGslMatrixClass::operator*=()",
                    "scaling failed");
  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator/=(double a)
{
  this->resetLU();
  *this *= (1./a);

  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator+=(const uqGslMatrixClass& rhs)
{
  this->resetLU();
  int iRC;
  iRC = gsl_matrix_add(m_mat,rhs.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
                    "uqGslMatrixClass::operator+=()",
                    "failed");

  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator-=(const uqGslMatrixClass& rhs)
{
  this->resetLU();
  int iRC;
  iRC = gsl_matrix_sub(m_mat,rhs.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
                    "uqGslMatrixClass::operator-=()",
                    "failed");

  return *this;
}

double&
uqGslMatrixClass::operator()(unsigned int i, unsigned int j)
{
  this->resetLU();
  if ((i >= m_mat->size1) ||
      (j >= m_mat->size2)) {
    std::cerr << "In uqGslMatrixClass::operator(i,j)"
              << ": i = " << i
              << ", j = " << j
              << ", m_mat->size1 = " << m_mat->size1
              << ", m_mat->size2 = " << m_mat->size2
              << std::endl;
    UQ_FATAL_TEST_MACRO(i >= m_mat->size1,
                        m_env.worldRank(),
                        "uqGslMatrixClass::operator(i,j)",
                        "i is too large");
    UQ_FATAL_TEST_MACRO(j >= m_mat->size2,
                        m_env.worldRank(),
                        "uqGslMatrixClass::operator(i,j)",
                        "j is too large");
  }
  return *gsl_matrix_ptr(m_mat,i,j);
}

const double&
uqGslMatrixClass::operator()(unsigned int i, unsigned int j) const
{
  UQ_FATAL_TEST_MACRO(i >= m_mat->size1,
                      m_env.worldRank(),
                      "uqGslMatrixClass::operator(i,j) const",
                      "i is too large");
  UQ_FATAL_TEST_MACRO(j >= m_mat->size2,
                      m_env.worldRank(),
                      "uqGslMatrixClass::operator(i,j) const",
                      "j is too large");
  return *gsl_matrix_const_ptr(m_mat,i,j);
}

void
uqGslMatrixClass::copy(const uqGslMatrixClass& src)
{
  this->resetLU();
  int iRC;
  iRC = gsl_matrix_memcpy(this->m_mat, src.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
                    "uqGslMatrixClass::copy()",
                    "failed");

  return;
}

void
uqGslMatrixClass::resetLU()
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  if (m_inverse) {
    delete m_inverse;
    m_inverse = NULL;
  }
  if (m_svdColMap) {
    delete m_svdColMap;
  }
  if (m_svdUmat) {
    delete m_svdUmat;
  }
  if (m_svdSvec) {
    delete m_svdSvec;
  }
  if (m_svdVmat) {
    delete m_svdVmat;
  }
  if (m_svdVTmat) {
    delete m_svdVTmat;
  }
  m_determinant   = -INFINITY;
  m_lnDeterminant = -INFINITY;
  if (m_permutation) {
    gsl_permutation_free(m_permutation);
    m_permutation = NULL;
  }
  m_signum = 0;

  return;
}

unsigned int
uqGslMatrixClass::numRowsLocal() const
{
  return m_mat->size1;
}

unsigned int
uqGslMatrixClass::numRowsGlobal() const
{
  return m_mat->size1;
}

unsigned int
uqGslMatrixClass::numCols() const
{
  return m_mat->size2;
}

double
uqGslMatrixClass::normFrob() const
{
  double value = 0.;

  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();
  double aux = 0.;
  for (unsigned int i = 0; i < nRows; i++) {
    for (unsigned int j = 0; j < nCols; j++) {
      aux = (*this)(i,j);
      value += aux*aux;
    }
  }

  return sqrt(value);
}

double
uqGslMatrixClass::normMax() const
{
  double value = 0.;

  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();
  double aux = 0.;
  for (unsigned int i = 0; i < nRows; i++) {
    for (unsigned int j = 0; j < nCols; j++) {
      aux = fabs((*this)(i,j));
      if (aux > value) value = aux;
    }
  }

  return value;
}

double
uqGslMatrixClass::max() const
{
  double value = -INFINITY;

  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();
  double aux = 0.;
  for (unsigned int i = 0; i < nRows; i++) {
    for (unsigned int j = 0; j < nCols; j++) {
      aux = (*this)(i,j);
      if (aux > value) value = aux;
    }
  }

  return value;
}

void
uqGslMatrixClass::cwSet(double value)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  for (unsigned int row = 0; row < nRows; ++row) {
    for (unsigned int col = 0; col < nCols; ++col) {
      *gsl_matrix_ptr(m_mat,row,col) = value;
    }
  }

  return;
}

void
uqGslMatrixClass::cwSet(
  unsigned int rowId,
  unsigned int colId,
  const uqGslMatrixClass& mat)
{
  UQ_FATAL_TEST_MACRO(rowId >= this->numRowsLocal(),
                      m_env.worldRank(),
                      "uqGslMatrixClass::cwSet()",
                      "invalid rowId");

  UQ_FATAL_TEST_MACRO((rowId + mat.numRowsLocal()) > this->numRowsLocal(),
                      m_env.worldRank(),
                      "uqGslMatrixClass::cwSet()",
                      "invalid vec.numRowsLocal()");

  UQ_FATAL_TEST_MACRO(colId >= this->numCols(),
                      m_env.worldRank(),
                      "uqGslMatrixClass::cwSet()",
                      "invalid colId");

  UQ_FATAL_TEST_MACRO((colId + mat.numCols()) > this->numCols(),
                      m_env.worldRank(),
                      "uqGslMatrixClass::cwSet()",
                      "invalid vec.numCols()");

  for (unsigned int i = 0; i < mat.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < mat.numCols(); ++j) {
      (*this)(rowId+i,colId+j) = mat(i,j);
    }
  }

  return;
}

int
uqGslMatrixClass::chol()
{
  int iRC;
  //std::cout << "Calling gsl_linalg_cholesky_decomp()..." << std::endl;
  gsl_error_handler_t* oldHandler;
  oldHandler = gsl_set_error_handler_off();
  iRC = gsl_linalg_cholesky_decomp(m_mat);
  if (iRC != 0) {
    std::cout << "In uqGslMatrixClass::chol()"
              << ": iRC = " << iRC
              << ", gsl error message = " << gsl_strerror(iRC)
              << std::endl;
  }
  gsl_set_error_handler(oldHandler);
  //std::cout << "Returned from gsl_linalg_cholesky_decomp() with iRC = " << iRC << std::endl;
  UQ_RC_MACRO(iRC, // Yes, *not* a fatal check on RC
              m_env.worldRank(),
              "uqGslMatrixClass::chol()",
              "matrix is not positive definite",
              UQ_MATRIX_IS_NOT_POS_DEFINITE_RC);

  return iRC;
}

int
uqGslMatrixClass::svd(uqGslMatrixClass& matU, uqGslVectorClass& vecS, uqGslMatrixClass& matVt) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((matU.numRowsLocal() != nRows) || (matU.numCols() != nCols),
                      m_env.worldRank(),
                      "uqGslMatrixClass::svd()",
                      "invalid matU");

  UQ_FATAL_TEST_MACRO((vecS.sizeLocal() != nCols), //std::min(nRows,nRows)),
                      m_env.worldRank(),
                      "uqGslMatrixClass::svd()",
                      "invalid vecS");

  UQ_FATAL_TEST_MACRO((matVt.numRowsLocal() != nCols) || (matVt.numCols() != nCols),
                      m_env.worldRank(),
                      "uqGslMatrixClass::svd()",
                      "invalid matVt");

  int iRC = internalSvd();

  matU  = *m_svdUmat;
  vecS  = *m_svdSvec;
  matVt = *m_svdVTmat;

  return iRC;
}

int
uqGslMatrixClass::svdSolve(const uqGslVectorClass& rhsVec, uqGslVectorClass& solVec) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((rhsVec.sizeLocal() != nRows),
                      m_env.worldRank(),
                      "uqGslMatrixClass::svdSolve()",
                      "invalid rhsVec");

  UQ_FATAL_TEST_MACRO((solVec.sizeLocal() != nCols),
                      m_env.worldRank(),
                      "uqGslMatrixClass::svdSolve()",
                      "invalid solVec");

  int iRC = internalSvd();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In uqGslMatrixClass::svdSolve():"
                            << "\n this->numRowsLocal()      = " << this->numRowsLocal()
                            << ", this->numCols()      = "       << this->numCols()
                            << "\n m_svdUmat->numRowsLocal() = " << m_svdUmat->numRowsLocal()
                            << ", m_svdUmat->numCols() = "       << m_svdUmat->numCols()
                            << "\n m_svdVmat->numRowsLocal() = " << m_svdVmat->numRowsLocal()
                            << ", m_svdVmat->numCols() = "       << m_svdVmat->numCols()
                            << "\n m_svdSvec->sizeLocal()    = " << m_svdSvec->sizeLocal()
                            << "\n rhsVec.sizeLocal()        = " << rhsVec.sizeLocal()
                            << "\n solVec.sizeLocal()        = " << solVec.sizeLocal()
                            << std::endl;
  }

  if (iRC == 0) iRC = gsl_linalg_SV_solve(m_svdUmat->data(), m_svdVmat->data(), m_svdSvec->data(), rhsVec.data(), solVec.data());

  return iRC;
}

int
uqGslMatrixClass::svdSolve(const uqGslMatrixClass& rhsMat, uqGslMatrixClass& solMat) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((rhsMat.numRowsLocal() != nRows),
                      m_env.worldRank(),
                      "uqGslMatrixClass::svdSolve()",
                      "invalid rhsMat");

  UQ_FATAL_TEST_MACRO((solMat.numRowsLocal() != nCols),
                      m_env.worldRank(),
                      "uqGslMatrixClass::svdSolve()",
                      "invalid solMat");

  UQ_FATAL_TEST_MACRO((rhsMat.numCols() != solMat.numCols()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::svdSolve()",
                      "rhsMat and solMat are not compatible");

  uqGslVectorClass rhsVec(m_env,rhsMat.map());
  uqGslVectorClass solVec(m_env,solMat.map());
  int iRC = 0;
  for (unsigned int j = 0; j < rhsMat.numCols(); ++j) {
    rhsVec = rhsMat.getColumn(j);
    iRC = this->svdSolve(rhsVec, solVec);
    if (iRC) break;
    solMat.setColumn(j,solVec);
  }

  return iRC;
}

const uqGslMatrixClass&
uqGslMatrixClass::svdMatU() const
{
  int iRC = 0;
  iRC = internalSvd();

  return *m_svdUmat;
}

const uqGslMatrixClass&
uqGslMatrixClass::svdMatV() const
{
  int iRC = 0;
  iRC = internalSvd();

  return *m_svdVmat;
}

int
uqGslMatrixClass::internalSvd() const
{
  int iRC = 0;

  if (m_svdColMap == NULL) {
    unsigned int nRows = this->numRowsLocal();
    unsigned int nCols = this->numCols();
    UQ_FATAL_TEST_MACRO(nRows < nCols,
                        m_env.worldRank(),
                        "uqGslMatrixClass::internalSvd()",
                        "GSL only supports cases where nRows >= nCols");

    m_svdColMap = new uqMapClass(this->numCols(),0,this->map().Comm()); // see 'uqVectorSpaceClass<.,.>::newMap()' in src/basic/src/uqGslVectorSpace.C
    m_svdUmat   = new uqGslMatrixClass(*this); // Yes, 'this'
    m_svdSvec   = new uqGslVectorClass(m_env,*m_svdColMap);
    m_svdVmat   = new uqGslMatrixClass(*m_svdSvec);
    m_svdVTmat  = new uqGslMatrixClass(*m_svdSvec);
    
  //uqGslVectorClass vecWork(*m_svdSvec );
    std::cout << "In uqGslMatrixClass::internalSvd()"
              << ", calling gsl_linalg_SV_decomp_jacobi()..."
              << ": nRows = " << nRows
              << ", nCols = " << nCols
              << std::endl;
    struct timeval timevalBegin;
    gettimeofday(&timevalBegin, NULL);
    gsl_error_handler_t* oldHandler;
    oldHandler = gsl_set_error_handler_off();
#if 1
    iRC = gsl_linalg_SV_decomp_jacobi(m_svdUmat->data(), m_svdVmat->data(), m_svdSvec->data());
#else
    iRC = gsl_linalg_SV_decomp(m_svdUmat->data(), m_svdVmat->data(), m_svdSvec->data(), vecWork.data());
#endif
    if (iRC != 0) {
      std::cout << "In uqGslMatrixClass::internalSvd()"
                << ": iRC = " << iRC
                << ", gsl error message = " << gsl_strerror(iRC)
                << std::endl;
    }
    gsl_set_error_handler(oldHandler);

    struct timeval timevalNow;
    gettimeofday(&timevalNow, NULL);
    std::cout << "In uqGslMatrixClass::internalSvd()"
              << ": returned from gsl_linalg_SV_decomp_jacobi() with iRC = " << iRC
              << " after " << timevalNow.tv_sec - timevalBegin.tv_sec
              << " seconds"
              << std::endl;
    UQ_RC_MACRO(iRC, // Yes, *not* a fatal check on RC
                m_env.worldRank(),
                "uqGslMatrixClass::internalSvd()",
                "matrix svd failed",
                UQ_MATRIX_SVD_FAILED_RC);
    *m_svdVTmat = m_svdVmat->transpose();
  }

  return iRC;
}

void
uqGslMatrixClass::zeroLower(bool includeDiagonal)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((nRows != nCols),
                      m_env.worldRank(),
                      "uqGslMatrixClass::zeroLower()",
                      "routine works only for square matrices");

  this->resetLU();

  if (includeDiagonal) {
    for (unsigned int i = 0; i < nRows; i++) {
      for (unsigned int j = 0; j <= i; j++) {
        (*this)(i,j) = 0.;
      }
    }
  }
  else {
    for (unsigned int i = 0; i < nRows; i++) {
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
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((nRows != nCols),
                      m_env.worldRank(),
                      "uqGslMatrixClass::zeroUpper()",
                      "routine works only for square matrices");

  this->resetLU();

  if (includeDiagonal) {
    for (unsigned int i = 0; i < nRows; i++) {
      for (unsigned int j = i; j < nCols; j++) {
        (*this)(i,j) = 0.;
      }
    }
  }
  else {
    for (unsigned int i = 0; i < nRows; i++) {
      for (unsigned int j = (i+1); j < nCols; j++) {
        (*this)(i,j) = 0.;
      }
    }
  }

  return;
}

void
uqGslMatrixClass::filterSmallValues(double thresholdValue)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();
  for (unsigned int i = 0; i < nRows; ++i) {
    for (unsigned int j = 0; j < nCols; ++j) {
      double aux = (*this)(i,j);
      // If 'thresholdValue' is negative, no values will be filtered
      if ((aux             < 0. ) &&
          (-thresholdValue < aux)) {
        (*this)(i,j) = 0.;
      }      
      if ((aux            > 0. ) &&
          (thresholdValue > aux)) {
        (*this)(i,j) = 0.;
      }      
    }
  }

  return;
}

void
uqGslMatrixClass::filterLargeValues(double thresholdValue)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();
  for (unsigned int i = 0; i < nRows; ++i) {
    for (unsigned int j = 0; j < nCols; ++j) {
      double aux = (*this)(i,j);
      // If 'thresholdValue' is negative, no values will be filtered
      if ((aux             < 0. ) &&
          (-thresholdValue > aux)) {
        (*this)(i,j) = 0.;
      }      
      if ((aux            > 0. ) &&
          (thresholdValue < aux)) {
        (*this)(i,j) = 0.;
      }      
    }
  }

  return;
}

uqGslMatrixClass
uqGslMatrixClass::transpose() const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((nRows != nCols),
                      m_env.worldRank(),
                      "uqGslMatrixClass::transpose()",
                      "routine works only for square matrices");

  uqGslMatrixClass mat(m_env,m_map,nCols);
  for (unsigned int row = 0; row < nRows; ++row) {
    for (unsigned int col = 0; col < nCols; ++col) {
      mat(row,col) = (*this)(col,row);
    }
  }

  return mat;
}

uqGslMatrixClass
uqGslMatrixClass::inverse() const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((nRows != nCols),
                      m_env.worldRank(),
                      "uqGslMatrixClass::inverse()",
                      "matrix is not square");

  if (m_inverse == NULL) {
    m_inverse = new uqGslMatrixClass(m_env,m_map,nCols);
    uqGslVectorClass unitVector(m_env,m_map);
    unitVector.cwSet(0.);
    uqGslVectorClass multVector(m_env,m_map);
    for (unsigned int j = 0; j < nCols; ++j) {
      if (j > 0) unitVector[j-1] = 0.;
      unitVector[j] = 1.;
      this->invertMultiply(unitVector, multVector);
      for (unsigned int i = 0; i < nRows; ++i) {
        (*m_inverse)(i,j) = multVector[i];
      }
    }
  }
  if (m_env.checkingLevel() >= 1) {
    *m_env.subDisplayFile() << "CHECKING In uqGslMatrixClass::inverse()"
                            << ": M.lnDet = "      << this->lnDeterminant()
                            << ", M^{-1}.lnDet = " << m_inverse->lnDeterminant()
                            << std::endl;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In uqGslMatrixClass::inverse():"
                            << "\n M = "        << *this
                            << "\n M^{-1} = "   << *m_inverse
                            << "\n M*M^{-1} = " << (*this)*(*m_inverse)
                            << "\n M^{-1}*M = " << (*m_inverse)*(*this)
                            << std::endl;
  }

  return *m_inverse;
}

void
uqGslMatrixClass::fillWithBlocksDiagonally(const std::vector<const uqGslMatrixClass* >& matrices)
{
  unsigned int sumNumRowsLocals = 0;
  unsigned int sumNumCols       = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    sumNumRowsLocals += matrices[i]->numRowsLocal();
    sumNumCols       += matrices[i]->numCols();
  }
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != sumNumRowsLocals,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithBlocksDiagonally(const)",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != sumNumCols,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithBlocksDiagonally(const)",
                      "inconsistent number of cols");

  unsigned int cumulativeRowId = 0;
  unsigned int cumulativeColId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(cumulativeRowId + rowId, cumulativeColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeRowId += nRows;
    cumulativeColId += nCols;
  }

  return;
}

void
uqGslMatrixClass::fillWithBlocksDiagonally(const std::vector<uqGslMatrixClass* >& matrices)
{
  unsigned int sumNumRowsLocals = 0;
  unsigned int sumNumCols       = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    sumNumRowsLocals += matrices[i]->numRowsLocal();
    sumNumCols       += matrices[i]->numCols();
  }
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != sumNumRowsLocals,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithBlocksDiagonally()",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != sumNumCols,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithBlocksDiagonally()",
                      "inconsistent number of cols");

  unsigned int cumulativeRowId = 0;
  unsigned int cumulativeColId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(cumulativeRowId + rowId, cumulativeColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    } 
    cumulativeRowId += nRows;
    cumulativeColId += nCols;
  }

  return;
}

void
uqGslMatrixClass::fillWithBlocksHorizontally(const std::vector<const uqGslMatrixClass* >& matrices)
{
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    UQ_FATAL_TEST_MACRO(this->numRowsLocal() != matrices[i]->numRowsLocal(),
                        m_env.worldRank(),
                        "uqGslMatrixClass::fillWithBlocksHorizontally(const)",
                        "inconsistent local number of rows");
    sumNumCols += matrices[i]->numCols();
  }
  UQ_FATAL_TEST_MACRO(this->numCols() != sumNumCols,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithBlocksHorizontally(const)",
                      "inconsistent number of cols");

  unsigned int cumulativeColId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(rowId, cumulativeColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeColId += nCols;
  }

  return;
}

void
uqGslMatrixClass::fillWithBlocksHorizontally(const std::vector<uqGslMatrixClass* >& matrices)
{
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    UQ_FATAL_TEST_MACRO(this->numRowsLocal() != matrices[i]->numRowsLocal(),
                        m_env.worldRank(),
                        "uqGslMatrixClass::fillWithBlocksHorizontally()",
                        "inconsistent local number of rows");
    sumNumCols += matrices[i]->numCols();
  }
  UQ_FATAL_TEST_MACRO(this->numCols() != sumNumCols,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithBlocksHorizontally()",
                      "inconsistent number of cols");

  unsigned int cumulativeColId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(rowId, cumulativeColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeColId += nCols;
  }

  return;
}

void
uqGslMatrixClass::fillWithBlocksVertically(const std::vector<const uqGslMatrixClass* >& matrices)
{
  unsigned int sumNumRows = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    UQ_FATAL_TEST_MACRO(this->numCols() != matrices[i]->numCols(),
                        m_env.worldRank(),
                        "uqGslMatrixClass::fillWithBlocksVertically(const)",
                        "inconsistent local number of cols");
    sumNumRows += matrices[i]->numRowsLocal();
  }
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != sumNumRows,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithBlocksVertically(const)",
                      "inconsistent number of rows");

  unsigned int cumulativeRowId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(cumulativeRowId + rowId, colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeRowId += nRows;
  }

  return;
}

void
uqGslMatrixClass::fillWithBlocksVertically(const std::vector<uqGslMatrixClass* >& matrices)
{
  unsigned int sumNumRows = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    UQ_FATAL_TEST_MACRO(this->numCols() != matrices[i]->numCols(),
                        m_env.worldRank(),
                        "uqGslMatrixClass::fillWithBlocksVertically()",
                        "inconsistent local number of cols");
    sumNumRows += matrices[i]->numRowsLocal();
  }
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != sumNumRows,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithBlocksVertically()",
                      "inconsistent number of rows");

  unsigned int cumulativeRowId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(cumulativeRowId + rowId, colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeRowId += nRows;
  }

  return;
}

void
uqGslMatrixClass::fillWithTensorProduct(const uqGslMatrixClass& mat1, const uqGslMatrixClass& mat2)
{
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != (mat1.numRowsLocal() * mat2.numRowsLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillTensorProduct(mat and mat)",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != (mat1.numCols() * mat2.numCols()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillTensorProduct(mat and mat)",
                      "inconsistent number of columns");

  for (unsigned int rowId1 = 0; rowId1 < mat1.numRowsLocal(); ++rowId1) {
    for (unsigned int colId1 = 0; colId1 < mat1.numCols(); ++colId1) {
      double multiplicativeFactor = mat1(rowId1,colId1);
      unsigned int targetRowId = rowId1 * mat2.numRowsLocal();
      unsigned int targetColId = colId1 * mat2.numCols();
      for (unsigned int rowId2 = 0; rowId2 < mat2.numRowsLocal(); ++rowId2) {
        for (unsigned int colId2 = 0; colId2 < mat2.numCols(); ++colId2) {
          (*this)(targetRowId + rowId2, targetColId + colId2) = multiplicativeFactor * mat2(rowId2,colId2);
        }
      }
    }
  } 

  return;
}

void
uqGslMatrixClass::fillWithTensorProduct(const uqGslMatrixClass& mat1, const uqGslVectorClass& vec2)
{
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != (mat1.numRowsLocal() * vec2.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillTensorProduct(mat and vec)",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != (mat1.numCols() * 1),
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillTensorProduct(mat and vec)",
                      "inconsistent number of columns");

  for (unsigned int rowId1 = 0; rowId1 < mat1.numRowsLocal(); ++rowId1) {
    for (unsigned int colId1 = 0; colId1 < mat1.numCols(); ++colId1) {
      double multiplicativeFactor = mat1(rowId1,colId1);
      unsigned int targetRowId = rowId1 * vec2.sizeLocal();
      unsigned int targetColId = colId1 * 1;
      for (unsigned int rowId2 = 0; rowId2 < vec2.sizeLocal(); ++rowId2) {
        for (unsigned int colId2 = 0; colId2 < 1; ++colId2) {
          (*this)(targetRowId + rowId2, targetColId + colId2) = multiplicativeFactor * vec2[rowId2];
        }
      }
    }
  } 


  return;
}

void
uqGslMatrixClass::fillWithTranspose(const uqGslMatrixClass& mat)
{
  unsigned int nRows = mat.numRowsLocal();
  unsigned int nCols = mat.numCols();
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != nCols,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithTranspose()",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != nRows,
                      m_env.worldRank(),
                      "uqGslMatrixClass::fillWithTranspose()",
                      "inconsistent number of cols");

  for (unsigned int row = 0; row < nRows; ++row) {
    for (unsigned int col = 0; col < nCols; ++col) {
      (*this)(col,row) = mat(row,col);
    }
  }

  return;
}

double
uqGslMatrixClass::determinant() const
{
  if (m_determinant == -INFINITY) {
    if (m_LU == NULL) {
      uqGslVectorClass tmpB(m_env,m_map);
      uqGslVectorClass tmpX(m_env,m_map);
      this->invertMultiply(tmpB,tmpX);
    }
    m_determinant   = gsl_linalg_LU_det(m_LU,m_signum); 
    m_lnDeterminant = gsl_linalg_LU_lndet(m_LU);
  }

  return m_determinant;
}

double
uqGslMatrixClass::lnDeterminant() const
{
  if (m_lnDeterminant == -INFINITY) {
    if (m_LU == NULL) {
      uqGslVectorClass tmpB(m_env,m_map);
      uqGslVectorClass tmpX(m_env,m_map);
      this->invertMultiply(tmpB,tmpX);
    }
    m_determinant   = gsl_linalg_LU_det(m_LU,m_signum); 
    m_lnDeterminant = gsl_linalg_LU_lndet(m_LU);
  }

  return m_lnDeterminant;
}

unsigned int
uqGslMatrixClass::rank(double absoluteZeroThreshold, double relativeZeroThreshold) const
{
  int iRC = 0;
  iRC = internalSvd();

  uqGslVectorClass relativeVec(*m_svdSvec);
  if (relativeVec[0] > 0.) {
    relativeVec = (1./relativeVec[0])*relativeVec;
  }

  unsigned int rankValue = 0;
  for (unsigned int i = 0; i < relativeVec.sizeLocal(); ++i) {
    if (( (*m_svdSvec)[i] >= absoluteZeroThreshold ) &&
        ( relativeVec [i] >= relativeZeroThreshold )) {
       rankValue += 1;
    }
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In uqGslMatrixClass::rank()"
                            << ": this->numRowsLocal() = "  << this->numRowsLocal()
                            << ", this->numCols() = "       << this->numCols()
                            << ", absoluteZeroThreshold = " << absoluteZeroThreshold
                            << ", relativeZeroThreshold = " << relativeZeroThreshold
                            << ", rankValue = "             << rankValue
                            << ", m_svdSvec = "             << *m_svdSvec
                            << ", relativeVec = "           << relativeVec
                            << std::endl;
  }

  return rankValue;
}

uqGslVectorClass
uqGslMatrixClass::multiply(
  const uqGslVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != x.sizeLocal()),
                      m_env.worldRank(),
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
  UQ_FATAL_TEST_MACRO((this->numCols() != x.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::multiply(), vector return void",
                      "matrix and x have incompatible sizes");

  UQ_FATAL_TEST_MACRO((this->numRowsLocal() != y.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::multiply(), vector return void",
                      "matrix and y have incompatible sizes");

  unsigned int sizeX = this->numCols();
  unsigned int sizeY = this->numRowsLocal();
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
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.worldRank(),
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
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::invertMultiply(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.sizeLocal() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::invertMultiply(), return void",
                      "solution and rhs have incompatible sizes");

  int iRC;
  if (m_LU == NULL) {
    UQ_FATAL_TEST_MACRO((m_permutation != NULL),
                        m_env.worldRank(),
                        "uqGslMatrixClass::invertMultiply()",
                        "m_permutation should be NULL");

    m_LU = gsl_matrix_calloc(this->numRowsLocal(),this->numCols());
    UQ_FATAL_TEST_MACRO((m_LU == NULL),
                        m_env.worldRank(),
                        "uqGslMatrixClass::invertMultiply()",
                        "gsl_matrix_calloc() failed");

    iRC = gsl_matrix_memcpy(m_LU, m_mat);
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.worldRank(),
                      "uqGslMatrixClass::invertMultiply()",
                      "gsl_matrix_memcpy() failed");

    m_permutation = gsl_permutation_calloc(numCols());
    UQ_FATAL_TEST_MACRO((m_permutation == NULL),
                        m_env.worldRank(),
                        "uqGslMatrixClass::invertMultiply()",
                        "gsl_permutation_calloc() failed");

    if (m_inDebugMode) {
      std::cout << "In uqGslMatrixClass::invertMultiply()"
                << ": before LU decomposition, m_LU = ";
      gsl_matrix_fprintf(stdout, m_LU, "%f");
      std::cout << std::endl;
    }

    iRC = gsl_linalg_LU_decomp(m_LU,m_permutation,&m_signum); 
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.worldRank(),
                      "uqGslMatrixClass::invertMultiply()",
                      "gsl_linalg_LU_decomp() failed");

    if (m_inDebugMode) {
      std::cout << "In uqGslMatrixClass::invertMultiply()"
                << ": after LU decomposition, m_LU = ";
      gsl_matrix_fprintf(stdout, m_LU, "%f");
      std::cout << std::endl;
    }
  }

  iRC = gsl_linalg_LU_solve(m_LU,m_permutation,b.data(),x.data()); 
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
                    "uqGslMatrixClass::invertMultiply()",
                    "gsl_linalg_LU_solve() failed");

  if (m_inDebugMode) {
    uqGslVectorClass tmpVec(b - (*this)*x);
    std::cout << "In uqGslMatrixClass::invertMultiply()"
              << ": ||b - Ax||_2 = "         << tmpVec.norm2()
              << ": ||b - Ax||_2/||b||_2 = " << tmpVec.norm2()/b.norm2()
              << std::endl;
  }

  return;
}

uqGslMatrixClass
uqGslMatrixClass::invertMultiply(const uqGslMatrixClass& B) const
{
  uqGslMatrixClass X(m_env,m_map,B.numCols());
  this->invertMultiply(B,X);

  return X;
}

void              
uqGslMatrixClass::invertMultiply(const uqGslMatrixClass& B, uqGslMatrixClass& X) const
{
  
  // Sanity Checks
  UQ_FATAL_RC_MACRO(((B.numRowsLocal() != X.numRowsLocal()) ||
		     (B.numCols()      != X.numCols()     )),
                    m_env.worldRank(),
		    "uqGslMatrixClass::invertMultiply()",
		    "Matrices B and X are incompatible");

  
  UQ_FATAL_RC_MACRO((this->numRowsLocal() != X.numRowsLocal()),
                    m_env.worldRank(),
		    "uqGslMatrixClass::invertMultiply()",
		    "This and X matrices are incompatible");

  // Some local variables used within the loop.
  uqGslVectorClass b(m_env, m_map);
  uqGslVectorClass x(m_env, m_map);

  for( unsigned int j = 0; j < B.numCols(); ++j )
    {
      b = B.getColumn( j );

      //invertMultiply will only do the LU once and store it. So we don't
      //need to worry about it doing LU multiple times.
      x = this->invertMultiply( b );

      X.setColumn( j, x );
    }

  return;
}

uqGslVectorClass
uqGslMatrixClass::invertMultiplyForceLU(
  const uqGslVectorClass& b) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::invertMultiply(), return vector",
                      "matrix and rhs have incompatible sizes");

  uqGslVectorClass x(m_env,m_map);
  this->invertMultiplyForceLU(b,x);

  return x;
}

void
uqGslMatrixClass::invertMultiplyForceLU(
  const uqGslVectorClass& b,
        uqGslVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::invertMultiplyForceLU(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.sizeLocal() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslMatrixClass::invertMultiplyForceLU(), return void",
                      "solution and rhs have incompatible sizes");

  int iRC;

  if ( m_LU == NULL ) {
    UQ_FATAL_TEST_MACRO((m_permutation != NULL),
                        m_env.worldRank(),
                        "uqGslMatrixClass::invertMultiplyForceLU()",
                        "m_permutation should be NULL");
    m_LU = gsl_matrix_calloc(this->numRowsLocal(),this->numCols());
  }
  UQ_FATAL_TEST_MACRO((m_LU == NULL),
		      m_env.worldRank(),
		      "uqGslMatrixClass::invertMultiplyForceLU()",
		      "gsl_matrix_calloc() failed");
  
  iRC = gsl_matrix_memcpy(m_LU, m_mat);
  UQ_FATAL_RC_MACRO(iRC,
		    m_env.worldRank(),
		    "uqGslMatrixClass::invertMultiplyForceLU()",
		    "gsl_matrix_memcpy() failed");
  
  if( m_permutation == NULL ) m_permutation = gsl_permutation_calloc(numCols());
  UQ_FATAL_TEST_MACRO((m_permutation == NULL),
		      m_env.worldRank(),
		      "uqGslMatrixClass::invertMultiplyForceLU()",
		      "gsl_permutation_calloc() failed");
  
  iRC = gsl_linalg_LU_decomp(m_LU,m_permutation,&m_signum); 
  UQ_FATAL_RC_MACRO(iRC,
		    m_env.worldRank(),
		    "uqGslMatrixClass::invertMultiplyForceLU()",
		    "gsl_linalg_LU_decomp() failed");

  iRC = gsl_linalg_LU_solve(m_LU,m_permutation,b.data(),x.data()); 
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
                    "uqGslMatrixClass::invertMultiplyForceLU()",
                    "gsl_linalg_LU_solve() failed");

  return;
}

void
uqGslMatrixClass::eigen(uqGslVectorClass& eigenValues, uqGslMatrixClass* eigenVectors) const
{
  unsigned int n = eigenValues.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      env().fullRank(),
                      "uqGslMatrixClass::eigen()",
                      "invalid input vector size");

  if (eigenVectors) {
    UQ_FATAL_TEST_MACRO((eigenValues.sizeLocal() != eigenVectors->numRowsLocal()),
                        env().fullRank(),
                        "uqGslVectorClass::eigen()",
                        "different input vector sizes");
  }

  if (eigenVectors == NULL) {
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc((size_t) n);
    gsl_eigen_symm(m_mat,eigenValues.data(),w);
    gsl_eigen_symm_free(w);
  }
  else {
    gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc((size_t) n);
    gsl_eigen_symmv(m_mat,eigenValues.data(),eigenVectors->m_mat,w);
    gsl_eigen_symmv_sort(eigenValues.data(),eigenVectors->m_mat,GSL_EIGEN_SORT_VAL_ASC);
    gsl_eigen_symmv_free(w);
  }

  return;
}

void
uqGslMatrixClass::largestEigen(double& eigenValue, uqGslVectorClass& eigenVector) const
{

  // Sanity Checks
  unsigned int n = eigenVector.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      env().fullRank(),
                      "uqGslMatrixClass::largestEigen()",
                      "invalid input vector size");

  /* The following notation is used:
     z = vector used in iteration that ends up being the eigenvector corresponding to the
         largest eigenvalue
     w = vector used in iteration that we extract the largest eigenvalue from.  */

  // Some parameters associated with the algorithm
  // TODO: Do we want to add the ability to have these set by the user?
  const unsigned int max_num_iterations = 10000;
  const double tolerance = 1.0e-13;

  // Create temporary working vectors.
  // TODO: We shouldn't have to use these - we ought to be able to work directly
  // TODO: with eigenValue and eigenVector. Optimize this once we have regression
  // TODO: tests.
  uqGslVectorClass z(m_env, m_map, 1.0 ); // Needs to be initialized to 1.0
  uqGslVectorClass w(m_env, m_map);

  // Some variables we use inside the loop.
  int index;
  double residual;
  double lambda;

  for( unsigned int k = 0; k < max_num_iterations; ++k )
    {
      w = (*this) * z;

      // For this algorithm, it's crucial to get the maximum in
      // absolute value, but then to normalize by the actual value
      // and *not* the absolute value.
      index = (w.abs()).getMaxValueIndex();

      lambda = w[index];

      z = (1.0/lambda) * w;

      // Here we use the norm of the residual as our convergence check:
      // norm( A*x - \lambda*x )
      residual = ( (*this)*z - lambda*z ).norm2();
      
      if( residual < tolerance )
	{
	  eigenValue = lambda;

	  // TODO: Do we want to normalize this so eigenVector.norm2() = 1?
	  eigenVector = z;
	  return;
	}
	
    }

  // If we reach this point, then we didn't converge. Print error message
  // to this effect.
  // TODO: We know we failed at this point - need a way to not need a test
  // TODO: Just specify the error.
  UQ_FATAL_TEST_MACRO((residual >= tolerance),
                      env().fullRank(),
                      "uqGslMatrixClass::largestEigen()",
                      "Maximum num iterations exceeded");

  
  return;
}

void
uqGslMatrixClass::smallestEigen(double& eigenValue, uqGslVectorClass& eigenVector) const
{
  // Sanity Checks
  unsigned int n = eigenVector.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      env().fullRank(),
                      "uqGslMatrixClass::smallestEigen()",
                      "invalid input vector size");

  /* The following notation is used:
     z = vector used in iteration that ends up being the eigenvector corresponding to the
         largest eigenvalue
     w = vector used in iteration that we extract the largest eigenvalue from.  */

  // Some parameters associated with the algorithm
  // TODO: Do we want to add the ability to have these set by the user?
  const unsigned int max_num_iterations = 1000;
  const double tolerance = 1.0e-13;

  // Create temporary working vectors.
  // TODO: We shouldn't have to use these - we ought to be able to work directly
  // TODO: with eigenValue and eigenVector. Optimize this once we have regression
  // TODO: tests.
  uqGslVectorClass z(m_env, m_map, 1.0 ); // Needs to be initialized to 1.0
  uqGslVectorClass w(m_env, m_map);

  // Some variables we use inside the loop.
  int index;
  double residual;
  double one_over_lambda;
  double lambda;

  for( unsigned int k = 0; k < max_num_iterations; ++k )
    {
      w = (*this).invertMultiplyForceLU( z );

      // For this algorithm, it's crucial to get the maximum in
      // absolute value, but then to normalize by the actual value
      // and *not* the absolute value.
      index = (w.abs()).getMaxValueIndex();

      // Remember: Inverse power method solves for max 1/lambda ==> lambda smallest
      one_over_lambda = w[index];

      lambda = 1.0/one_over_lambda;
      
      z = lambda * w;

      // Here we use the norm of the residual as our convergence check:
      // norm( A*x - \lambda*x )
      residual = ( (*this)*z - lambda*z ).norm2();
      
      if( residual < tolerance )
	{
	  eigenValue = lambda;

	  // TODO: Do we want to normalize this so eigenVector.norm2() = 1?
	  eigenVector = z;
	  return;
	}
	
    }

  // If we reach this point, then we didn't converge. Print error message
  // to this effect.
  // TODO: We know we failed at this point - need a way to not need a test
  // TODO: Just specify the error.
  UQ_FATAL_TEST_MACRO((residual >= tolerance),
                      env().fullRank(),
                      "uqGslMatrixClass::smallestEigen()",
                      "Maximum num iterations exceeded");

  return;
}

void
uqGslMatrixClass::getColumn(unsigned int column_num, uqGslVectorClass& column) const
{
  // Sanity checks
  UQ_FATAL_TEST_MACRO(column_num >= this->numCols(),
                      env().fullRank(),
                      "uqGslMatrixClass::getColumn",
                      "Specified row number not within range");

  UQ_FATAL_TEST_MACRO((column.sizeLocal() != this->numRowsLocal()),
                      env().fullRank(),
                      "uqGslMatrixClass::getColumn",
                      "column vector not same size as this matrix");

  // Temporary working vector
  gsl_vector* gsl_column = gsl_vector_alloc( column.sizeLocal() );

  int error_code = gsl_matrix_get_col( gsl_column, m_mat, column_num );
  UQ_FATAL_RC_MACRO( (error_code != 0),
		     env().fullRank(),
		     "uqGslMatrixClass::getColumn()",
		     "gsl_matrix_get_col failed");

  // Copy column from gsl matrix into our GslVector object
  for( unsigned int i = 0; i < column.sizeLocal(); ++i )
    {
      column[i] = gsl_vector_get( gsl_column, i );
    }

  // Clean up gsl temporaries
  gsl_vector_free( gsl_column );

  return;
}

void
uqGslMatrixClass::getRow(unsigned int row_num, uqGslVectorClass& row) const
{
  // Sanity checks
  UQ_FATAL_TEST_MACRO(row_num >= this->numRowsLocal(),
                      env().fullRank(),
                      "uqGslMatrixClass::getRow",
                      "Specified row number not within range");

  UQ_FATAL_TEST_MACRO((row.sizeLocal() != this->numCols()),
                      env().fullRank(),
                      "uqGslMatrixClass::getRow",
                      "row vector not same size as this matrix");

  // Temporary working vector
  gsl_vector* gsl_row = gsl_vector_alloc( row.sizeLocal() );

  int error_code = gsl_matrix_get_row( gsl_row, m_mat, row_num );
  UQ_FATAL_RC_MACRO( (error_code != 0),
		     env().fullRank(),
		     "uqGslMatrixClass::getRow()",
		     "gsl_matrix_get_row failed");

  // Copy row from gsl matrix into our GslVector object
  for( unsigned int i = 0; i < row.sizeLocal(); ++i )
    {
      row[i] = gsl_vector_get( gsl_row, i );
    }

  // Clean up gsl temporaries
  gsl_vector_free( gsl_row );

  return;
}

uqGslVectorClass
uqGslMatrixClass::getRow(unsigned int row_num ) const
{
  // We rely on the void getRow's to do sanity checks
  // So we don't do them here.

  uqGslVectorClass row(m_env, m_map);

  this->getRow( row_num, row );
  
  return row;
}

uqGslVectorClass
uqGslMatrixClass::getColumn(unsigned int column_num ) const
{
  // We rely on the void getRow's to do sanity checks
  // So we don't do them here.

  uqGslVectorClass column(m_env, m_map);

  this->getColumn( column_num, column );
  
  return column;
}

void
uqGslMatrixClass::setRow(unsigned int row_num, const uqGslVectorClass& row)
{
  this->resetLU();
  // Sanity checks
  UQ_FATAL_TEST_MACRO(row_num >= this->numRowsLocal(),
                      env().fullRank(),
                      "uqGslMatrixClass::setRow",
                      "Specified row number not within range");

  UQ_FATAL_TEST_MACRO((row.sizeLocal() != this->numCols()),
                      env().fullRank(),
                      "uqGslMatrixClass::setRow",
                      "row vector not same size as this matrix");

  gsl_vector* gsl_row = row.data();

  int error_code = gsl_matrix_set_row( m_mat, row_num, gsl_row );
  UQ_FATAL_RC_MACRO( (error_code != 0),
		     env().fullRank(),
		     "uqGslMatrixClass::setRow()",
		     "gsl_matrix_set_row failed");

  return;
}

void
uqGslMatrixClass::setColumn(unsigned int column_num, const uqGslVectorClass& column)
{
  this->resetLU();
  // Sanity checks
  UQ_FATAL_TEST_MACRO(column_num >= this->numCols(),
                      env().fullRank(),
                      "uqGslMatrixClass::setColumn",
                      "Specified column number not within range");

  UQ_FATAL_TEST_MACRO((column.sizeLocal() != this->numRowsLocal()),
                      env().fullRank(),
                      "uqGslMatrixClass::setColumn",
                      "column vector not same size as this matrix");

  gsl_vector* gsl_column = column.data();

  int error_code = gsl_matrix_set_col( m_mat, column_num, gsl_column );
  UQ_FATAL_RC_MACRO( (error_code != 0),
		     env().fullRank(),
		     "uqGslMatrixClass::setColumn()",
		     "gsl_matrix_set_col failed");

  return;
}

void
uqGslMatrixClass::mpiSum( const uqMpiCommClass& comm, uqGslMatrixClass& M_global ) const
{
  // Sanity Checks
  UQ_FATAL_RC_MACRO(((this->numRowsLocal() != M_global.numRowsLocal()) ||
                     (this->numCols()      != M_global.numCols()     )),
		    env().fullRank(),
		    "uqGslMatrixClass::mpiSum()",
		    "local and global matrices incompatible");

  /* TODO: Probably a better way to handle this unpacking/packing of data */
  int size = M_global.numRowsLocal()*M_global.numCols();
  std::vector<double> local( size, 0.0 );
  std::vector<double> global( size, 0.0 );

  int k;
  for( unsigned int i = 0; i < this->numRowsLocal(); i++ )
    {
      for( unsigned int j = 0; j < this->numCols(); j++ )
	{
	  k = i + j*M_global.numCols();
	  
	  local[k] = (*this)(i,j);
	}
    }

  comm.Allreduce((void*) &local[0], (void*) &global[0], size, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                 "uqGslMatrixClass::mpiSum()",
                 "failed MPI.Allreduce()");

  for( unsigned int i = 0; i < this->numRowsLocal(); i++ )
    {
      for( unsigned int j = 0; j < this->numCols(); j++ )
	{
	  k = i + j*M_global.numCols();
	  
	  M_global(i,j) = global[k];
	}
    }

  return;
}

void
uqGslMatrixClass::matlabLinearInterpExtrap(
  const uqGslVectorClass& x1Vec,
  const uqGslMatrixClass& y1Mat,
  const uqGslVectorClass& x2Vec)
{
  UQ_FATAL_TEST_MACRO(x1Vec.sizeLocal() <= 1,
                      m_env.worldRank(),
                      "uqGslMatrixClass::matlabLinearInterpExtrap()",
                      "invalid 'x1' size");

  UQ_FATAL_TEST_MACRO(x1Vec.sizeLocal() != y1Mat.numRowsLocal(),
                      m_env.worldRank(),
                      "uqGslMatrixClass::matlabLinearInterpExtrap()",
                      "invalid 'x1' and 'y1' sizes");

  UQ_FATAL_TEST_MACRO(x2Vec.sizeLocal() != this->numRowsLocal(),
                      m_env.worldRank(),
                      "uqGslMatrixClass::matlabLinearInterpExtrap()",
                      "invalid 'x2' and 'this' sizes");

  UQ_FATAL_TEST_MACRO(y1Mat.numCols() != this->numCols(),
                      m_env.worldRank(),
                      "uqGslMatrixClass::matlabLinearInterpExtrap()",
                      "invalid 'y1' and 'this' sizes");

  uqGslVectorClass y1Vec(x1Vec);
  uqGslVectorClass y2Vec(x2Vec);
  for (unsigned int colId = 0; colId < y1Mat.numCols(); ++colId) {
    y1Mat.getColumn(colId,y1Vec);
    y2Vec.matlabLinearInterpExtrap(x1Vec,y1Vec,x2Vec);
    this->setColumn(colId,y2Vec);
  }

  return;
}

gsl_matrix*
uqGslMatrixClass::data()
{
  return m_mat;
}

void
uqGslMatrixClass::print(std::ostream& os) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  if (m_printHorizontally) {
    for (unsigned int i = 0; i < nRows; ++i) {
      for (unsigned int j = 0; j < nCols; ++j) {
        os << (*this)(i,j)
           << " ";
      }
      if (i != (nRows-1)) os << "; ";
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

void
uqGslMatrixClass::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqGslMatrixClass::subWriteContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "uqGslMatrixClass::subWriteContents()",
                      "implemented just for sequential vectors for now");

  uqFilePtrSetStruct filePtrSet;
  if (m_env.openOutputFile(fileName,
                           fileType, // "m or hdf"
                           allowedSubEnvIds,
                           false,
                           filePtrSet)) {
    unsigned int nRows = this->numRowsLocal();
    unsigned int nCols = this->numCols();
    *filePtrSet.ofsVar << varNamePrefix << "_sub" << m_env.subIdString() << " = zeros(" << nRows
                       << ","                                                           << nCols
                       << ");"
                       << std::endl;
    *filePtrSet.ofsVar << varNamePrefix << "_sub" << m_env.subIdString() << " = [";

    for (unsigned int i = 0; i < nRows; ++i) {
      for (unsigned int j = 0; j < nCols; ++j) {
        *filePtrSet.ofsVar << (*this)(i,j)
                << " ";
      }
      *filePtrSet.ofsVar << "\n";
    }
    *filePtrSet.ofsVar << "];\n";

    m_env.closeFile(filePtrSet,fileType);
  }

  return;
}

void
uqGslMatrixClass::subReadContents(
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds)
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqGslMatrixClass::subReadContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "uqGslMatrixClass::subReadContents()",
                      "implemented just for sequential vectors for now");

  uqFilePtrSetStruct filePtrSet;
  if (m_env.openInputFile(fileName,
                          fileType, // "m or hdf"
                          allowedSubEnvIds,
                          filePtrSet)) {

     // palms
    double nRowsLocal = this->numRowsLocal();

    // In the logic below, the id of a line' begins with value 0 (zero)
    unsigned int idOfMyFirstLine = 1;
    unsigned int idOfMyLastLine = nRowsLocal;
    unsigned int nCols = this->numCols();

    // Read number of matrix rows in the file by taking care of the first line,
    // which resembles something like 'variable_name = zeros(n_rows,n_cols);'
    std::string tmpString;

    // Read 'variable name' string
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Just read '" << tmpString << "'" << std::endl;

    // Read '=' sign
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Just read '" << tmpString << "'" << std::endl;
    UQ_FATAL_TEST_MACRO(tmpString != "=",
                        m_env.worldRank(),
                        "uqGslMatrixClass::subReadContents()",
                        "string should be the '=' sign");

    // Read 'zeros(n_rows,n_cols)' string
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Just read '" << tmpString << "'" << std::endl;
    unsigned int posInTmpString = 6;

    // Isolate 'n_rows' in a string
    char nRowsString[tmpString.size()-posInTmpString+1];
    unsigned int posInRowsString = 0;
    do {
      UQ_FATAL_TEST_MACRO(posInTmpString >= tmpString.size(),
                          m_env.worldRank(),
                          "uqGslMatrixClass::subReadContents()",
                          "symbol ',' not found in first line of file");
      nRowsString[posInRowsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ',');
    nRowsString[posInRowsString] = '\0';

    // Isolate 'n_cols' in a string
    posInTmpString++; // Avoid reading ',' char
    char nColsString[tmpString.size()-posInTmpString+1];
    unsigned int posInColsString = 0;
    do {
      UQ_FATAL_TEST_MACRO(posInTmpString >= tmpString.size(),
                          m_env.worldRank(),
                          "uqGslMatrixClass::subReadContents()",
                          "symbol ')' not found in first line of file");
      nColsString[posInColsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ')');
    nColsString[posInColsString] = '\0';

    // Convert 'n_rows' and 'n_cols' strings to numbers
    unsigned int numRowsInFile = (unsigned int) strtod(nRowsString,NULL);
    unsigned int numColsInFile = (unsigned int) strtod(nColsString,NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGslMatrixClass::subReadContents()"
                              << ": fullRank "        << m_env.fullRank()
                              << ", numRowsInFile = " << numRowsInFile
                              << ", numColsInFile = " << numColsInFile
                              << ", nRowsLocal = "    << nRowsLocal
                              << ", nCols = "         << nCols
                              << std::endl;
    }

    // Check if [num of rows in file] == [requested matrix row size]
    UQ_FATAL_TEST_MACRO(numRowsInFile != nRowsLocal,
                        m_env.worldRank(),
                        "uqGslMatrixClass::subReadContents()",
                        "size of vec in file is not big enough");

    // Check if [num of cols in file] == [num cols in current matrix]
    UQ_FATAL_TEST_MACRO(numColsInFile != nCols,
                        m_env.worldRank(),
                        "uqGslMatrixClass::subReadContents()",
                        "number of parameters of vec in file is different than number of parameters in this vec object");

    // Code common to any core in a communicator
    unsigned int maxCharsPerLine = 64*nCols; // Up to about 60 characters to represent each parameter value

    unsigned int lineId = 0;
    while (lineId < idOfMyFirstLine) {
      filePtrSet.ifsVar->ignore(maxCharsPerLine,'\n');
      lineId++;
    };

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGslMatrixClass::subReadContents()"
                              << ": beginning to read input actual data"
                              << std::endl;
    }

    // Take care of initial part of the first data line,
    // which resembles something like 'variable_name = [value1 value2 ...'
    // Read 'variable name' string
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Core 0 just read '" << tmpString << "'" << std::endl;

    // Read '=' sign
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Core 0 just read '" << tmpString << "'" << std::endl;
    UQ_FATAL_TEST_MACRO(tmpString != "=",
                        m_env.worldRank(),
                        "uqGslMatrixClass::subReadContents()",
                        "in core 0, string should be the '=' sign");

    // Take into account the ' [' portion
    std::streampos tmpPos = filePtrSet.ifsVar->tellg();
    filePtrSet.ifsVar->seekg(tmpPos+(std::streampos)2);

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGslMatrixClass::subReadContents()"
                              << ": beginning to read lines with numbers only"
                              << ", lineId = "          << lineId
                              << ", idOfMyFirstLine = " << idOfMyFirstLine
                              << ", idOfMyLastLine = "  << idOfMyLastLine
                              << std::endl;
    }

    double tmpRead;
    while (lineId <= idOfMyLastLine) {
      for (unsigned int j = 0; j < nCols; ++j) {
        *filePtrSet.ifsVar >> tmpRead;
        (*this)(lineId-idOfMyFirstLine,j) = tmpRead;
      }
      lineId++;
    };

    m_env.closeFile(filePtrSet,fileType);
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
  unsigned int m1Rows = m1.numRowsLocal();
  unsigned int m1Cols = m1.numCols();
  unsigned int m2Rows = m2.numRowsLocal();
  unsigned int m2Cols = m2.numCols();

  UQ_FATAL_TEST_MACRO((m1Cols != m2Rows),
                      m1.env().worldRank(),
                      "uqGslMatrixClass operator*(matrix,matrix)",
                      "different sizes m1Cols and m2Rows");

  uqGslMatrixClass mat(m1.env(),m1.map(),m2Cols);

  //std::cout << "In uqGslMatrixClass(mat * mat): m1Cols = " << m1Cols << std::endl;

  unsigned int commonSize = m1Cols;
  for (unsigned int row1 = 0; row1 < m1Rows; ++row1) {
    for (unsigned int col2 = 0; col2 < m2Cols; ++col2) {
      double result = 0.;
      for (unsigned int k = 0; k < commonSize; ++k) {
        result += m1(row1,k)*m2(k,col2);
      }
      mat(row1,col2) = result;
    }
    //std::cout << "In uqGslMatrixClass(mat * mat): ended row " << row1 << std::endl;
  }

  return mat;
}

uqGslMatrixClass operator+(const uqGslMatrixClass& m1, const uqGslMatrixClass& m2)
{
  uqGslMatrixClass answer(m1);
  answer += m2;
  return answer;
}

uqGslMatrixClass operator-(const uqGslMatrixClass& m1, const uqGslMatrixClass& m2)
{
  uqGslMatrixClass answer(m1);
  answer -= m2;
  return answer;
}

uqGslMatrixClass matrixProduct(const uqGslVectorClass& v1, const uqGslVectorClass& v2)
{
  unsigned int nRows = v1.sizeLocal();
  unsigned int nCols = v2.sizeLocal();
  uqGslMatrixClass answer(v1.env(),v1.map(),nCols);

  for (unsigned int i = 0; i < nRows; ++i) {
    double value1 = v1[i];
    for (unsigned int j = 0; j < nCols; ++j) {
      answer(i,j) = value1*v2[j];
    }
  }

  return answer;
}

uqGslMatrixClass leftDiagScaling(const uqGslVectorClass& vec, const uqGslMatrixClass& mat)
{
  unsigned int vSize = vec.sizeLocal();
  unsigned int mRows = mat.numRowsLocal();
  unsigned int mCols = mat.numCols();

  UQ_FATAL_TEST_MACRO((vSize != mRows),
                      mat.env().worldRank(),
                      "uqGslMatrixClass leftDiagScaling(vector,matrix)",
                      "size of vector is different from the number of rows in matrix");

  UQ_FATAL_TEST_MACRO((mCols != mRows),
                      mat.env().worldRank(),
                      "uqGslMatrixClass leftDiagScaling(vector,matrix)",
                      "routine currently works for square matrices only");

  uqGslMatrixClass answer(mat);
  for (unsigned int i = 0; i < mRows; ++i) {
    double vecValue = vec[i];
    for (unsigned int j = 0; j < mCols; ++j) {
      answer(i,j) *= vecValue;
    }
  }

  return answer;
}

uqGslMatrixClass rightDiagScaling(const uqGslMatrixClass& mat, const uqGslVectorClass& vec)
{
  unsigned int vSize = vec.sizeLocal();
  unsigned int mRows = mat.numRowsLocal();
  unsigned int mCols = mat.numCols();

  UQ_FATAL_TEST_MACRO((vSize != mCols),
                      mat.env().worldRank(),
                      "uqGslMatrixClass rightDiagScaling(matrix,vector)",
                      "size of vector is different from the number of cols in matrix");

  UQ_FATAL_TEST_MACRO((mCols != mRows),
                      mat.env().worldRank(),
                      "uqGslMatrixClass rightDiagScaling(matrix,vector)",
                      "routine currently works for square matrices only");

  uqGslMatrixClass answer(mat);
  for (unsigned int j = 0; j < mCols; ++j) {
    double vecValue = vec[j];
    for (unsigned int i = 0; i < mRows; ++i) {
      answer(i,j) *= vecValue;
    }
  }

  return answer;
}
