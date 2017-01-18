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

#include <queso/GslMatrix.h>
#include <queso/GslVector.h>
#include <queso/Defines.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <sys/time.h>
#include <cmath>

namespace QUESO {

GslMatrix::GslMatrix( // can be a rectangular matrix
  const BaseEnvironment& env,
  const Map&             map,
  unsigned int                  nCols)
  :
  Matrix  (env,map),
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
  m_signum       (0),
  m_isSingular   (false)
{
  queso_require_msg(m_mat, "null matrix generated");
}

// Shaped constructor --------------------------------------------------
GslMatrix::GslMatrix( // square matrix
  const BaseEnvironment& env,
  const Map&             map,
  double                        diagValue)
  :
  Matrix  (env,map),
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
  m_signum       (0),
  m_isSingular   (false)
{
  queso_require_msg(m_mat, "null matrix generated");

  for (unsigned int i = 0; i < m_mat->size1; ++i) {
    (*this)(i,i) = diagValue;
  }
}
// Shaped constructor --------------------------------------------------
GslMatrix::GslMatrix( // square matrix
  const GslVector& v,
  double                  diagValue)
  :
  Matrix  (v.env(),v.map()),
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
  m_signum       (0),
  m_isSingular   (false)
{
  queso_require_msg(m_mat, "null matrix generated");

  for (unsigned int i = 0; i < m_mat->size1; ++i) {
    (*this)(i,i) = diagValue;
  }
}
// Shaped constructor --------------------------------------------------
GslMatrix::GslMatrix(const GslVector& v) // square matrix
  :
  Matrix  (v.env(),v.map()),
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
  m_signum       (0),
  m_isSingular   (false)
{
  queso_require_msg(m_mat, "null matrix generated");

  unsigned int dim = v.sizeLocal();
  for (unsigned int i = 0; i < dim; ++i) {
    (*this)(i,i) = v[i];
  }
}

// Shaped constructor --------------------------------------------------
GslMatrix::GslMatrix(const GslMatrix& B) // can be a rectangular matrix
  :
  Matrix  (B.env(),B.map()),
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
  m_signum       (0),
  m_isSingular   (false)
{
  queso_require_msg(m_mat, "null vector generated");
  this->Matrix::base_copy(B);
  this->copy(B);
}

// Destructor ----------------------------------------------------------
GslMatrix::~GslMatrix()
{
  this->resetLU();
  if (m_mat) gsl_matrix_free(m_mat);
}

// Set methods (operators) ---------------------------------------------
GslMatrix&
GslMatrix::operator=(const GslMatrix& obj)
{
  this->resetLU();
  this->copy(obj);
  return *this;
}

GslMatrix&
GslMatrix::operator*=(double a)
{
  this->resetLU();
  int iRC;
  iRC = gsl_matrix_scale(m_mat,a);
  queso_require_msg(!(iRC), "scaling failed");
  return *this;
}

GslMatrix&
GslMatrix::operator/=(double a)
{
  this->resetLU();
  *this *= (1./a);

  return *this;
}

GslMatrix&
GslMatrix::operator+=(const GslMatrix& rhs)
{
  this->resetLU();
  int iRC;
  iRC = gsl_matrix_add(m_mat,rhs.m_mat);
  queso_require_msg(!(iRC), "failed");

  return *this;
}

GslMatrix&
GslMatrix::operator-=(const GslMatrix& rhs)
{
  this->resetLU();
  int iRC;
  iRC = gsl_matrix_sub(m_mat,rhs.m_mat);
  queso_require_msg(!(iRC), "failed");

  return *this;
}



void
GslMatrix::copy(const GslMatrix& src)
{
  // FIXME - should we be calling Matrix::base_copy here? - RHS
  this->resetLU();
  int iRC;
  iRC = gsl_matrix_memcpy(this->m_mat, src.m_mat);
  queso_require_msg(!(iRC), "failed");

  return;
}

void
GslMatrix::resetLU()
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
    m_svdColMap = NULL;
  }
  if (m_svdUmat) {
    delete m_svdUmat;
    m_svdUmat = NULL;
  }
  if (m_svdSvec) {
    delete m_svdSvec;
    m_svdSvec = NULL;
  }
  if (m_svdVmat) {
    delete m_svdVmat;
    m_svdVmat = NULL;
  }
  if (m_svdVTmat) {
    delete m_svdVTmat;
    m_svdVTmat = NULL;
  }
  m_determinant   = -INFINITY;
  m_lnDeterminant = -INFINITY;
  if (m_permutation) {
    gsl_permutation_free(m_permutation);
    m_permutation = NULL;
  }
  m_signum = 0;
  m_isSingular = false;

  return;
}

unsigned int
GslMatrix::numRowsLocal() const
{
  return m_mat->size1;
}

unsigned int
GslMatrix::numRowsGlobal() const
{
  return m_mat->size1;
}

unsigned int
GslMatrix::numCols() const
{
  return m_mat->size2;
}

double
GslMatrix::normFrob() const
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
GslMatrix::normMax() const
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
GslMatrix::max() const
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
GslMatrix::cwSet(double value)
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
GslMatrix::cwSet(
  unsigned int            initialTargetRowId,
  unsigned int            initialTargetColId,
  const GslMatrix& mat)
{
  queso_require_less_msg(initialTargetRowId, this->numRowsLocal(), "invalid initialTargetRowId");

  queso_require_less_equal_msg((initialTargetRowId + mat.numRowsLocal()), this->numRowsLocal(), "invalid vec.numRowsLocal()");

  queso_require_less_msg(initialTargetColId, this->numCols(), "invalid initialTargetColId");

  queso_require_less_equal_msg((initialTargetColId + mat.numCols()), this->numCols(), "invalid vec.numCols()");

  for (unsigned int i = 0; i < mat.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < mat.numCols(); ++j) {
      (*this)(initialTargetRowId+i,initialTargetColId+j) = mat(i,j);
    }
  }

  return;
}

void
GslMatrix::cwExtract(
  unsigned int      initialTargetRowId,
  unsigned int      initialTargetColId,
  GslMatrix& mat) const
{
  queso_require_less_msg(initialTargetRowId, this->numRowsLocal(), "invalid initialTargetRowId");

  queso_require_less_equal_msg((initialTargetRowId + mat.numRowsLocal()), this->numRowsLocal(), "invalid vec.numRowsLocal()");

  queso_require_less_msg(initialTargetColId, this->numCols(), "invalid initialTargetColId");

  queso_require_less_equal_msg((initialTargetColId + mat.numCols()), this->numCols(), "invalid vec.numCols()");

  for (unsigned int i = 0; i < mat.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < mat.numCols(); ++j) {
      mat(i,j) = (*this)(initialTargetRowId+i,initialTargetColId+j) ;
    }
  }

  return;
}

int
GslMatrix::chol()
{
  int iRC;
  //std::cout << "Calling gsl_linalg_cholesky_decomp()..." << std::endl;
  gsl_error_handler_t* oldHandler;
  oldHandler = gsl_set_error_handler_off();
  iRC = gsl_linalg_cholesky_decomp(m_mat);
  if (iRC != 0) {
    std::cerr << "In GslMatrix::chol()"
              << ": iRC = " << iRC
              << ", gsl error message = " << gsl_strerror(iRC)
              << std::endl;
    std::cerr << "Here is the offending matrix: " << std::endl;
    std::cerr << *this << std::endl;
  }
  gsl_set_error_handler(oldHandler);
  //std::cout << "Returned from gsl_linalg_cholesky_decomp() with iRC = " << iRC << std::endl;
  UQ_RC_MACRO(iRC, // Yes, *not* a fatal check on RC
              m_env.worldRank(),
              "GslMatrix::chol()",
              "matrix is not positive definite",
              UQ_MATRIX_IS_NOT_POS_DEFINITE_RC);

  return iRC;
}

int
GslMatrix::svd(GslMatrix& matU, GslVector& vecS, GslMatrix& matVt) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  queso_require_msg(!((matU.numRowsLocal() != nRows) || (matU.numCols() != nCols)), "invalid matU");

  queso_require_equal_to_msg(vecS.sizeLocal(), nCols, "invalid vecS");

  queso_require_msg(!((matVt.numRowsLocal() != nCols) || (matVt.numCols() != nCols)), "invalid matVt");

  int iRC = internalSvd();

  matU  = *m_svdUmat;
  vecS  = *m_svdSvec;
  matVt = *m_svdVTmat;

  return iRC;
}

int
GslMatrix::svdSolve(const GslVector& rhsVec, GslVector& solVec) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  queso_require_equal_to_msg(rhsVec.sizeLocal(), nRows, "invalid rhsVec");

  queso_require_equal_to_msg(solVec.sizeLocal(), nCols, "invalid solVec");

  int iRC = internalSvd();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In GslMatrix::svdSolve():"
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
GslMatrix::svdSolve(const GslMatrix& rhsMat, GslMatrix& solMat) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  queso_require_equal_to_msg(rhsMat.numRowsLocal(), nRows, "invalid rhsMat");

  queso_require_equal_to_msg(solMat.numRowsLocal(), nCols, "invalid solMat");

  queso_require_equal_to_msg(rhsMat.numCols(), solMat.numCols(), "rhsMat and solMat are not compatible");

  GslVector rhsVec(m_env,rhsMat.map());
  GslVector solVec(m_env,solMat.map());
  int iRC = 0;
  for (unsigned int j = 0; j < rhsMat.numCols(); ++j) {
    rhsVec = rhsMat.getColumn(j);
    iRC = this->svdSolve(rhsVec, solVec);
    if (iRC) break;
    solMat.setColumn(j,solVec);
  }

  return iRC;
}

const GslMatrix&
GslMatrix::svdMatU() const
{
  int iRC = 0;
  iRC = internalSvd();
  if (iRC) {}; // just to remove compiler warning

  return *m_svdUmat;
}

const GslMatrix&
GslMatrix::svdMatV() const
{
  int iRC = 0;
  iRC = internalSvd();
  if (iRC) {}; // just to remove compiler warning

  return *m_svdVmat;
}

int
GslMatrix::internalSvd() const
{
  int iRC = 0;

  if (m_svdColMap == NULL) {
    unsigned int nRows = this->numRowsLocal();
    unsigned int nCols = this->numCols();
    queso_require_greater_equal_msg(nRows, nCols, "GSL only supports cases where nRows >= nCols");

    m_svdColMap = new Map(this->numCols(),0,this->map().Comm()); // see 'VectorSpace<.,.>::newMap()' in src/basic/src/GslVectorSpace.C
    m_svdUmat   = new GslMatrix(*this); // Yes, 'this'
    m_svdSvec   = new GslVector(m_env,*m_svdColMap);
    m_svdVmat   = new GslMatrix(*m_svdSvec);
    m_svdVTmat  = new GslMatrix(*m_svdSvec);

    //std::cout << "In GslMatrix::internalSvd()"
    //          << ", calling gsl_linalg_SV_decomp_jacobi()..."
    //          << ": nRows = " << nRows
    //          << ", nCols = " << nCols
    //          << std::endl;
    struct timeval timevalBegin;
    gettimeofday(&timevalBegin, NULL);
    gsl_error_handler_t* oldHandler;
    oldHandler = gsl_set_error_handler_off();
#if 1
    iRC = gsl_linalg_SV_decomp_jacobi(m_svdUmat->data(), m_svdVmat->data(), m_svdSvec->data());
#else
    GslVector vecWork(*m_svdSvec );
    iRC = gsl_linalg_SV_decomp(m_svdUmat->data(), m_svdVmat->data(), m_svdSvec->data(), vecWork.data());
#endif
    if (iRC != 0) {
      std::cerr << "In GslMatrix::internalSvd()"
                << ": iRC = " << iRC
                << ", gsl error message = " << gsl_strerror(iRC)
                << std::endl;
    }
    gsl_set_error_handler(oldHandler);

    struct timeval timevalNow;
    gettimeofday(&timevalNow, NULL);
    //std::cout << "In GslMatrix::internalSvd()"
    //          << ": returned from gsl_linalg_SV_decomp_jacobi() with iRC = " << iRC
    //          << " after " << timevalNow.tv_sec - timevalBegin.tv_sec
    //          << " seconds"
    //          << std::endl;
    UQ_RC_MACRO(iRC, // Yes, *not* a fatal check on RC
                m_env.worldRank(),
                "GslMatrix::internalSvd()",
                "matrix svd failed",
                UQ_MATRIX_SVD_FAILED_RC);
    *m_svdVTmat = m_svdVmat->transpose();
  }

  return iRC;
}



void
GslMatrix::zeroLower(bool includeDiagonal)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  queso_require_equal_to_msg(nRows, nCols, "routine works only for square matrices");

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
GslMatrix::zeroUpper(bool includeDiagonal)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  queso_require_equal_to_msg(nRows, nCols, "routine works only for square matrices");

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
GslMatrix::filterSmallValues(double thresholdValue)
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
GslMatrix::filterLargeValues(double thresholdValue)
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

GslMatrix
GslMatrix::transpose() const
{
  unsigned int nRows = this->numRowsGlobal();
  unsigned int nCols = this->numCols();

  const MpiComm & comm = this->map().Comm();
  Map serial_map(nCols, 0, comm);

  GslMatrix mat(m_env,serial_map,nRows);

  for (unsigned int row = 0; row < nRows; ++row) {
    for (unsigned int col = 0; col < nCols; ++col) {
      mat(col,row) = (*this)(row,col);
    }
  }

  return mat;
}

GslMatrix
GslMatrix::inverse() const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  queso_require_equal_to_msg(nRows, nCols, "matrix is not square");

  if (m_inverse == NULL) {
    m_inverse = new GslMatrix(m_env,m_map,nCols);
    GslVector unitVector(m_env,m_map);
    unitVector.cwSet(0.);
    GslVector multVector(m_env,m_map);
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
    *m_env.subDisplayFile() << "CHECKING In GslMatrix::inverse()"
                            << ": M.lnDet = "      << this->lnDeterminant()
                            << ", M^{-1}.lnDet = " << m_inverse->lnDeterminant()
                            << std::endl;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In GslMatrix::inverse():"
                            << "\n M = "        << *this
                            << "\n M^{-1} = "   << *m_inverse
                            << "\n M*M^{-1} = " << (*this)*(*m_inverse)
                            << "\n M^{-1}*M = " << (*m_inverse)*(*this)
                            << std::endl;
  }

  return *m_inverse;
}

void
GslMatrix::fillWithBlocksDiagonally(
  unsigned int                                 initialTargetRowId,
  unsigned int                                 initialTargetColId,
  const std::vector<const GslMatrix* >& matrices,
  bool                                         checkForExactNumRowsMatching,
  bool                                         checkForExactNumColsMatching)
{
  unsigned int sumNumRowsLocals = 0;
  unsigned int sumNumCols       = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    sumNumRowsLocals += matrices[i]->numRowsLocal();
    sumNumCols       += matrices[i]->numCols();
  }
  queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + sumNumRowsLocals), "too big number of rows");
  if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + sumNumRowsLocals), "inconsistent number of rows");
  queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + sumNumCols), "too big number of cols");
  if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + sumNumCols), "inconsistent number of cols");

  unsigned int cumulativeRowId = 0;
  unsigned int cumulativeColId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(initialTargetRowId + cumulativeRowId + rowId, initialTargetColId + cumulativeColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeRowId += nRows;
    cumulativeColId += nCols;
  }

  return;
}

void
GslMatrix::fillWithBlocksDiagonally(
  unsigned int                           initialTargetRowId,
  unsigned int                           initialTargetColId,
  const std::vector<GslMatrix* >& matrices,
  bool                                   checkForExactNumRowsMatching,
  bool                                   checkForExactNumColsMatching)
{
  unsigned int sumNumRowsLocals = 0;
  unsigned int sumNumCols       = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    sumNumRowsLocals += matrices[i]->numRowsLocal();
    sumNumCols       += matrices[i]->numCols();
  }
  queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + sumNumRowsLocals), "too big number of rows");
  if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + sumNumRowsLocals), "inconsistent number of rows");
  queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + sumNumCols), "too big number of cols");
  if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + sumNumCols), "inconsistent number of cols");

  unsigned int cumulativeRowId = 0;
  unsigned int cumulativeColId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(initialTargetRowId + cumulativeRowId + rowId, initialTargetColId + cumulativeColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeRowId += nRows;
    cumulativeColId += nCols;
  }

  return;
}

void
GslMatrix::fillWithBlocksHorizontally(
  unsigned int                                 initialTargetRowId,
  unsigned int                                 initialTargetColId,
  const std::vector<const GslMatrix* >& matrices,
  bool                                         checkForExactNumRowsMatching,
  bool                                         checkForExactNumColsMatching)
{
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + matrices[i]->numRowsLocal()), "too big number of rows");
    if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + matrices[i]->numRowsLocal()), "inconsistent number of rows");
    sumNumCols += matrices[i]->numCols();
  }
  queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + sumNumCols), "too big number of cols");
  if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + sumNumCols), "inconsistent number of cols");

  unsigned int cumulativeColId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(initialTargetRowId + rowId, initialTargetColId + cumulativeColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeColId += nCols;
  }

  return;
}

void
GslMatrix::fillWithBlocksHorizontally(
  unsigned int                           initialTargetRowId,
  unsigned int                           initialTargetColId,
  const std::vector<GslMatrix* >& matrices,
  bool                                   checkForExactNumRowsMatching,
  bool                                   checkForExactNumColsMatching)
{
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + matrices[i]->numRowsLocal()), "too big number of rows");
    if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + matrices[i]->numRowsLocal()), "inconsistent number of rows");
    sumNumCols += matrices[i]->numCols();
  }
  queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + sumNumCols), "too big number of cols");
  if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + sumNumCols), "inconsistent number of cols");

  unsigned int cumulativeColId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(initialTargetRowId + rowId, initialTargetColId + cumulativeColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeColId += nCols;
  }

  return;
}

void
GslMatrix::fillWithBlocksVertically( // checar
  unsigned int                                 initialTargetRowId,
  unsigned int                                 initialTargetColId,
  const std::vector<const GslMatrix* >& matrices,
  bool                                         checkForExactNumRowsMatching,
  bool                                         checkForExactNumColsMatching)
{
  unsigned int sumNumRows = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + matrices[i]->numCols()), "too big number of cols");
    if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + matrices[i]->numCols()), "inconsistent number of cols");
    sumNumRows += matrices[i]->numRowsLocal();
  }
  queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + sumNumRows), "too big number of rows");
  if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + sumNumRows), "inconsistent number of rows");

  unsigned int cumulativeRowId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(initialTargetRowId + cumulativeRowId + rowId, initialTargetColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeRowId += nRows;
  }

  return;
}

void
GslMatrix::fillWithBlocksVertically( // checar
  unsigned int                           initialTargetRowId,
  unsigned int                           initialTargetColId,
  const std::vector<GslMatrix* >& matrices,
  bool                                   checkForExactNumRowsMatching,
  bool                                   checkForExactNumColsMatching)
{
  unsigned int sumNumRows = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + matrices[i]->numCols()), "too big number of cols");
    if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + matrices[i]->numCols()), "inconsistent number of cols");
    sumNumRows += matrices[i]->numRowsLocal();
  }
  queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + sumNumRows), "too big number of rows");
  if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + sumNumRows), "inconsistent number of rows");

  unsigned int cumulativeRowId = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    unsigned int nRows = matrices[i]->numRowsLocal();
    unsigned int nCols = matrices[i]->numCols();
    for (unsigned int rowId = 0; rowId < nRows; ++rowId) {
      for (unsigned int colId = 0; colId < nCols; ++colId) {
        (*this)(initialTargetRowId + cumulativeRowId + rowId, initialTargetColId + colId) = (*(matrices[i]))(rowId,colId);
      }
    }
    cumulativeRowId += nRows;
  }

  return;
}

void
GslMatrix::fillWithTensorProduct(
  unsigned int            initialTargetRowId,
  unsigned int            initialTargetColId,
  const GslMatrix& mat1,
  const GslMatrix& mat2,
  bool                    checkForExactNumRowsMatching,
  bool                    checkForExactNumColsMatching)
{
  queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + (mat1.numRowsLocal() * mat2.numRowsLocal())), "too big number of rows");
  if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + (mat1.numRowsLocal() * mat2.numRowsLocal())), "inconsistent number of rows");
  queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + (mat1.numCols() * mat2.numCols())), "too big number of columns");
  if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + (mat1.numCols() * mat2.numCols())), "inconsistent number of columns");

  for (unsigned int rowId1 = 0; rowId1 < mat1.numRowsLocal(); ++rowId1) {
    for (unsigned int colId1 = 0; colId1 < mat1.numCols(); ++colId1) {
      double multiplicativeFactor = mat1(rowId1,colId1);
      unsigned int targetRowId = rowId1 * mat2.numRowsLocal();
      unsigned int targetColId = colId1 * mat2.numCols();
      for (unsigned int rowId2 = 0; rowId2 < mat2.numRowsLocal(); ++rowId2) {
        for (unsigned int colId2 = 0; colId2 < mat2.numCols(); ++colId2) {
          (*this)(initialTargetRowId + targetRowId + rowId2, initialTargetColId + targetColId + colId2) = multiplicativeFactor * mat2(rowId2,colId2);
        }
      }
    }
  }

  return;
}

void
GslMatrix::fillWithTensorProduct(
  unsigned int            initialTargetRowId,
  unsigned int            initialTargetColId,
  const GslMatrix& mat1,
  const GslVector& vec2,
  bool                    checkForExactNumRowsMatching,
  bool                    checkForExactNumColsMatching)
{
  queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + (mat1.numRowsLocal() * vec2.sizeLocal())), "too big number of rows");
  if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + (mat1.numRowsLocal() * vec2.sizeLocal())), "inconsistent number of rows");
  queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + (mat1.numCols() * 1)), "too big number of columns");
  if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + (mat1.numCols() * 1)), "inconsistent number of columns");

  for (unsigned int rowId1 = 0; rowId1 < mat1.numRowsLocal(); ++rowId1) {
    for (unsigned int colId1 = 0; colId1 < mat1.numCols(); ++colId1) {
      double multiplicativeFactor = mat1(rowId1,colId1);
      unsigned int targetRowId = rowId1 * vec2.sizeLocal();
      unsigned int targetColId = colId1 * 1;
      for (unsigned int rowId2 = 0; rowId2 < vec2.sizeLocal(); ++rowId2) {
        for (unsigned int colId2 = 0; colId2 < 1; ++colId2) {
          (*this)(initialTargetRowId + targetRowId + rowId2, initialTargetColId + targetColId + colId2) = multiplicativeFactor * vec2[rowId2];
        }
      }
    }
  }


  return;
}

void
GslMatrix::fillWithTranspose(
  unsigned int            initialTargetRowId,
  unsigned int            initialTargetColId,
  const GslMatrix& mat,
  bool                    checkForExactNumRowsMatching,
  bool                    checkForExactNumColsMatching)
{
  unsigned int nRows = mat.numRowsLocal();
  unsigned int nCols = mat.numCols();
  queso_require_greater_equal_msg(this->numRowsLocal(), (initialTargetRowId + nCols), "too big number of rows");
  if (checkForExactNumRowsMatching) queso_require_equal_to_msg(this->numRowsLocal(), (initialTargetRowId + nCols), "inconsistent number of rows");
  queso_require_greater_equal_msg(this->numCols(), (initialTargetColId + nRows), "too big number of cols");
  if (checkForExactNumColsMatching) queso_require_equal_to_msg(this->numCols(), (initialTargetColId + nRows), "inconsistent number of cols");

  for (unsigned int row = 0; row < nRows; ++row) {
    for (unsigned int col = 0; col < nCols; ++col) {
      (*this)(initialTargetRowId + col, initialTargetColId + row) = mat(row,col);
    }
  }

  return;
}

double
GslMatrix::determinant() const
{
  if (m_determinant == -INFINITY) {
    if (m_LU == NULL) {
      GslVector tmpB(m_env,m_map);
      GslVector tmpX(m_env,m_map);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In GslMatrix::determinant()"
                                << ": before 'this->invertMultiply()'"
                                << std::endl;
      }
      this->invertMultiply(tmpB,tmpX);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In GslMatrix::determinant()"
                                << ": after 'this->invertMultiply()'"
                                << std::endl;
      }
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslMatrix::determinant()"
                              << ": before 'gsl_linalg_LU_det()'"
                              << std::endl;
    }
    m_determinant   = gsl_linalg_LU_det(m_LU,m_signum);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslMatrix::determinant()"
                              << ": after 'gsl_linalg_LU_det()'"
                              << std::endl;
    }
    m_lnDeterminant = gsl_linalg_LU_lndet(m_LU);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslMatrix::determinant()"
                              << ": after 'gsl_linalg_LU_lndet()'"
                              << std::endl;
    }
  }

  return m_determinant;
}

double
GslMatrix::lnDeterminant() const
{
  if (m_lnDeterminant == -INFINITY) {
    if (m_LU == NULL) {
      GslVector tmpB(m_env,m_map);
      GslVector tmpX(m_env,m_map);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In GslMatrix::lnDeterminant()"
                                << ": before 'this->invertMultiply()'"
                                << std::endl;
      }
      this->invertMultiply(tmpB,tmpX);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In GslMatrix::lnDeterminant()"
                                << ": after 'this->invertMultiply()'"
                                << std::endl;
      }
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslMatrix::lnDeterminant()"
                              << ": before 'gsl_linalg_LU_det()'"
                              << std::endl;
    }
    m_determinant   = gsl_linalg_LU_det(m_LU,m_signum);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslMatrix::lnDeterminant()"
                              << ": after 'gsl_linalg_LU_det()'"
                              << std::endl;
    }
    m_lnDeterminant = gsl_linalg_LU_lndet(m_LU);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslMatrix::lnDeterminant()"
                              << ": before 'gsl_linalg_LU_lndet()'"
                              << std::endl;
    }
  }

  return m_lnDeterminant;
}

unsigned int
GslMatrix::rank(double absoluteZeroThreshold, double relativeZeroThreshold) const
{
  int iRC = 0;
  iRC = internalSvd();
  if (iRC) {}; // just to remove compiler warning

  GslVector relativeVec(*m_svdSvec);
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
    *m_env.subDisplayFile() << "In GslMatrix::rank()"
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

GslVector
GslMatrix::multiply(
  const GslVector& x) const
{
  queso_require_equal_to_msg(this->numCols(), x.sizeLocal(), "matrix and vector have incompatible sizes");

  GslVector y(m_env,m_map);
  this->multiply(x,y);

  return y;
}

void
GslMatrix::multiply(
  const GslVector& x,
        GslVector& y) const
{
  queso_require_equal_to_msg(this->numCols(), x.sizeLocal(), "matrix and x have incompatible sizes");

  queso_require_equal_to_msg(this->numRowsLocal(), y.sizeLocal(), "matrix and y have incompatible sizes");

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

GslMatrix
GslMatrix::multiply(
  const GslMatrix & X) const
{
  GslMatrix Y(m_env,m_map,X.numCols());
  this->multiply(X,Y);

  return Y;
}



void
GslMatrix::multiply(
  const GslMatrix & X,
        GslMatrix & Y) const
{
  queso_require_equal_to_msg(this->numCols(), X.numRowsGlobal(), "matrix and X have incompatible sizes");
  queso_require_equal_to_msg(this->numRowsGlobal(), Y.numRowsGlobal(), "matrix and Y have incompatible sizes");
  queso_require_equal_to_msg(X.numCols(), Y.numCols(), "X and Y have incompatible sizes");

  const unsigned int m_s = this->numRowsGlobal();
  const unsigned int p_s = this->numCols();
  const unsigned int n_s = X.numCols();

  for (unsigned int k=0; k<p_s; k++)
    for (unsigned int j=0; j<n_s; j++)
      if (X(k,j) != 0.)
        for (unsigned int i=0; i<m_s; i++)
          Y(i,j) += (*this)(i,k) * X(k,j);
}


GslVector
GslMatrix::invertMultiply(
  const GslVector& b) const
{
  queso_require_equal_to_msg(this->numCols(), b.sizeLocal(), "matrix and rhs have incompatible sizes");

  GslVector x(m_env,m_map);
  this->invertMultiply(b,x);

  return x;
}

void
GslMatrix::invertMultiply(
  const GslVector& b,
        GslVector& x) const
{
  queso_require_equal_to_msg(this->numCols(), b.sizeLocal(), "matrix and rhs have incompatible sizes");

  queso_require_equal_to_msg(x.sizeLocal(), b.sizeLocal(), "solution and rhs have incompatible sizes");

  int iRC;
  if (m_LU == NULL) {
    queso_require_msg(!(m_permutation), "m_permutation should be NULL");

    m_LU = gsl_matrix_calloc(this->numRowsLocal(),this->numCols());
    queso_require_msg(m_LU, "gsl_matrix_calloc() failed");

    iRC = gsl_matrix_memcpy(m_LU, m_mat);
    queso_require_msg(!(iRC), "gsl_matrix_memcpy() failed");

    m_permutation = gsl_permutation_calloc(numCols());
    queso_require_msg(m_permutation, "gsl_permutation_calloc() failed");

    if (m_inDebugMode) {
      std::cout << "In GslMatrix::invertMultiply()"
                << ": before LU decomposition, m_LU = ";
      gsl_matrix_fprintf(stdout, m_LU, "%f");
      std::cout << std::endl;
    }

    gsl_error_handler_t* oldHandler;
    oldHandler = gsl_set_error_handler_off();
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslMatrix::invertMultiply()"
                              << ": before 'gsl_linalg_LU_decomp()'"
                              << std::endl;
    }
    iRC = gsl_linalg_LU_decomp(m_LU,m_permutation,&m_signum);
    if (iRC != 0) {
      std::cerr << "In GslMatrix::invertMultiply()"
                << ", after gsl_linalg_LU_decomp()"
                << ": iRC = " << iRC
                << ", gsl error message = " << gsl_strerror(iRC)
                << std::endl;
    }
    gsl_set_error_handler(oldHandler);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslMatrix::invertMultiply()"
                              << ": after 'gsl_linalg_LU_decomp()'"
                              << ", IRC = " << iRC
                              << std::endl;
    }
    queso_require_msg(!(iRC), "gsl_linalg_LU_decomp() failed");

    if (m_inDebugMode) {
      std::cout << "In GslMatrix::invertMultiply()"
                << ": after LU decomposition, m_LU = ";
      gsl_matrix_fprintf(stdout, m_LU, "%f");
      std::cout << std::endl;
    }
  }

  gsl_error_handler_t* oldHandler;
  oldHandler = gsl_set_error_handler_off();
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GslMatrix::invertMultiply()"
                            << ": before 'gsl_linalg_LU_solve()'"
                            << std::endl;
  }
  iRC = gsl_linalg_LU_solve(m_LU,m_permutation,b.data(),x.data());
  if (iRC != 0) {
    m_isSingular = true;
    std::cerr << "In GslMatrix::invertMultiply()"
              << ", after gsl_linalg_LU_solve()"
              << ": iRC = " << iRC
              << ", gsl error message = " << gsl_strerror(iRC)
              << std::endl;
  }
  gsl_set_error_handler(oldHandler);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In GslMatrix::invertMultiply()"
                            << ": after 'gsl_linalg_LU_solve()'"
                            << ", IRC = " << iRC
                            << std::endl;
  }


  //                  m_env.worldRank(),
  //                  "GslMatrix::invertMultiply()",
  //                  "gsl_linalg_LU_solve() failed");

  if (m_inDebugMode) {
    GslVector tmpVec(b - (*this)*x);
    std::cout << "In GslMatrix::invertMultiply()"
              << ": ||b - Ax||_2 = "         << tmpVec.norm2()
              << ": ||b - Ax||_2/||b||_2 = " << tmpVec.norm2()/b.norm2()
              << std::endl;
  }

  return;
}

GslMatrix
GslMatrix::invertMultiply(const GslMatrix& B) const
{
  GslMatrix X(m_env,m_map,B.numCols());
  this->invertMultiply(B,X);

  return X;
}

void
GslMatrix::invertMultiply(const GslMatrix& B, GslMatrix& X) const
{
  // Sanity Checks
  queso_require_equal_to_msg(B.numRowsLocal(), X.numRowsLocal(),
		             "Matrices B and X are incompatible");
  queso_require_equal_to_msg(B.numCols(),      X.numCols(),
		             "Matrices B and X are incompatible");
  queso_require_equal_to_msg(this->numRowsLocal(), X.numRowsLocal(),
                             "This and X matrices are incompatible");

  // Some local variables used within the loop.
  GslVector b(m_env, m_map);
  GslVector x(m_env, m_map);

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

GslVector
GslMatrix::invertMultiplyForceLU(
  const GslVector& b) const
{
  queso_require_equal_to_msg(this->numCols(), b.sizeLocal(), "matrix and rhs have incompatible sizes");

  GslVector x(m_env,m_map);
  this->invertMultiplyForceLU(b,x);

  return x;
}

void
GslMatrix::invertMultiplyForceLU(
  const GslVector& b,
        GslVector& x) const
{
  queso_require_equal_to_msg(this->numCols(), b.sizeLocal(), "matrix and rhs have incompatible sizes");

  queso_require_equal_to_msg(x.sizeLocal(), b.sizeLocal(), "solution and rhs have incompatible sizes");

  int iRC;

  if ( m_LU == NULL ) {
    queso_require_msg(!(m_permutation), "m_permutation should be NULL");
    m_LU = gsl_matrix_calloc(this->numRowsLocal(),this->numCols());
  }
  queso_require_msg(m_LU, "gsl_matrix_calloc() failed");

  iRC = gsl_matrix_memcpy(m_LU, m_mat);
  queso_require_msg(!(iRC), "gsl_matrix_memcpy() failed");

  if( m_permutation == NULL ) m_permutation = gsl_permutation_calloc(numCols());
  queso_require_msg(m_permutation, "gsl_permutation_calloc() failed");

  iRC = gsl_linalg_LU_decomp(m_LU,m_permutation,&m_signum);
  queso_require_msg(!(iRC), "gsl_linalg_LU_decomp() failed");

  iRC = gsl_linalg_LU_solve(m_LU,m_permutation,b.data(),x.data());
  if (iRC != 0) {
    m_isSingular = true;
  }
  queso_require_msg(!(iRC), "gsl_linalg_LU_solve() failed");

  return;
}

void
GslMatrix::eigen(GslVector& eigenValues, GslMatrix* eigenVectors) const
{
  unsigned int n = eigenValues.sizeLocal();

  queso_require_not_equal_to_msg(n, 0, "invalid input vector size");

  if (eigenVectors) {
    queso_require_equal_to_msg(eigenValues.sizeLocal(), eigenVectors->numRowsLocal(), "different input vector sizes");
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
GslMatrix::largestEigen(double& eigenValue, GslVector& eigenVector) const
{

  // Sanity Checks
  unsigned int n = eigenVector.sizeLocal();

  queso_require_not_equal_to_msg(n, 0, "invalid input vector size");

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
  GslVector z(m_env, m_map, 1.0 ); // Needs to be initialized to 1.0
  GslVector w(m_env, m_map);

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
  queso_require_less_msg(residual, tolerance, "Maximum num iterations exceeded");


  return;
}

void
GslMatrix::smallestEigen(double& eigenValue, GslVector& eigenVector) const
{
  // Sanity Checks
  unsigned int n = eigenVector.sizeLocal();

  queso_require_not_equal_to_msg(n, 0, "invalid input vector size");

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
  GslVector z(m_env, m_map, 1.0 ); // Needs to be initialized to 1.0
  GslVector w(m_env, m_map);

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
  queso_require_less_msg(residual, tolerance, "Maximum num iterations exceeded");

  return;
}

void
GslMatrix::getColumn(unsigned int column_num, GslVector& column) const
{
  // Sanity checks
  queso_require_less_msg(column_num, this->numCols(), "Specified row number not within range");

  queso_require_equal_to_msg(column.sizeLocal(), this->numRowsLocal(), "column vector not same size as this matrix");

  // Temporary working vector
  gsl_vector* gsl_column = gsl_vector_alloc( column.sizeLocal() );

  int error_code = gsl_matrix_get_col( gsl_column, m_mat, column_num );
  queso_require_equal_to_msg(error_code, 0, "gsl_matrix_get_col failed");

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
GslMatrix::getRow(unsigned int row_num, GslVector& row) const
{
  // Sanity checks
  queso_require_less_msg(row_num, this->numRowsLocal(), "Specified row number not within range");

  queso_require_equal_to_msg(row.sizeLocal(), this->numCols(), "row vector not same size as this matrix");

  // Temporary working vector
  gsl_vector* gsl_row = gsl_vector_alloc( row.sizeLocal() );

  int error_code = gsl_matrix_get_row( gsl_row, m_mat, row_num );
  queso_require_equal_to_msg(error_code, 0, "gsl_matrix_get_row failed");

  // Copy row from gsl matrix into our GslVector object
  for( unsigned int i = 0; i < row.sizeLocal(); ++i )
    {
      row[i] = gsl_vector_get( gsl_row, i );
    }

  // Clean up gsl temporaries
  gsl_vector_free( gsl_row );

  return;
}

GslVector
GslMatrix::getRow(unsigned int row_num ) const
{
  // We rely on the void getRow's to do sanity checks
  // So we don't do them here.

  GslVector row(m_env, m_map);

  this->getRow( row_num, row );

  return row;
}

GslVector
GslMatrix::getColumn(unsigned int column_num ) const
{
  // We rely on the void getRow's to do sanity checks
  // So we don't do them here.

  GslVector column(m_env, m_map);

  this->getColumn( column_num, column );

  return column;
}

void
GslMatrix::setRow(unsigned int row_num, const GslVector& row)
{
  this->resetLU();
  // Sanity checks
  queso_require_less_msg(row_num, this->numRowsLocal(), "Specified row number not within range");

  queso_require_equal_to_msg(row.sizeLocal(), this->numCols(), "row vector not same size as this matrix");

  gsl_vector* gsl_row = row.data();

  int error_code = gsl_matrix_set_row( m_mat, row_num, gsl_row );
  queso_require_equal_to_msg(error_code, 0, "gsl_matrix_set_row failed");

  return;
}

void
GslMatrix::setColumn(unsigned int column_num, const GslVector& column)
{
  this->resetLU();
  // Sanity checks
  queso_require_less_msg(column_num, this->numCols(), "Specified column number not within range");

  queso_require_equal_to_msg(column.sizeLocal(), this->numRowsLocal(), "column vector not same size as this matrix");

  gsl_vector* gsl_column = column.data();

  int error_code = gsl_matrix_set_col( m_mat, column_num, gsl_column );
  queso_require_equal_to_msg(error_code, 0, "gsl_matrix_set_col failed");

  return;
}

void
GslMatrix::mpiSum( const MpiComm& comm, GslMatrix& M_global ) const
{
  // Sanity Checks
  UQ_FATAL_RC_MACRO(((this->numRowsLocal() != M_global.numRowsLocal()) ||
                     (this->numCols()      != M_global.numCols()     )),
		    env().fullRank(),
		    "GslMatrix::mpiSum()",
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

  comm.Allreduce<double>(&local[0], &global[0], size, RawValue_MPI_SUM,
                 "GslMatrix::mpiSum()",
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
GslMatrix::matlabLinearInterpExtrap(
  const GslVector& x1Vec,
  const GslMatrix& y1Mat,
  const GslVector& x2Vec)
{
  queso_require_greater_msg(x1Vec.sizeLocal(), 1, "invalid 'x1' size");

  queso_require_equal_to_msg(x1Vec.sizeLocal(), y1Mat.numRowsLocal(), "invalid 'x1' and 'y1' sizes");

  queso_require_equal_to_msg(x2Vec.sizeLocal(), this->numRowsLocal(), "invalid 'x2' and 'this' sizes");

  queso_require_equal_to_msg(y1Mat.numCols(), this->numCols(), "invalid 'y1' and 'this' sizes");

  GslVector y1Vec(x1Vec);
  GslVector y2Vec(x2Vec);
  for (unsigned int colId = 0; colId < y1Mat.numCols(); ++colId) {
    y1Mat.getColumn(colId,y1Vec);
    y2Vec.matlabLinearInterpExtrap(x1Vec,y1Vec,x2Vec);
    this->setColumn(colId,y2Vec);
  }

  return;
}

gsl_matrix*
GslMatrix::data()
{
  return m_mat;
}

void
GslMatrix::print(std::ostream& os) const
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
GslMatrix::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  queso_require_greater_equal_msg(m_env.subRank(), 0, "unexpected subRank");

  queso_require_less_equal_msg(this->numOfProcsForStorage(), 1, "implemented just for sequential vectors for now");

  FilePtrSetStruct filePtrSet;
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
GslMatrix::subReadContents(
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds)
{
  queso_require_greater_equal_msg(m_env.subRank(), 0, "unexpected subRank");

  queso_require_less_equal_msg(this->numOfProcsForStorage(), 1, "implemented just for sequential vectors for now");

  FilePtrSetStruct filePtrSet;
  if (m_env.openInputFile(fileName,
                          fileType, // "m or hdf"
                          allowedSubEnvIds,
                          filePtrSet)) {

    // palms
    unsigned int nRowsLocal = this->numRowsLocal();

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
    queso_require_equal_to_msg(tmpString, std::string("="), std::string("string should be the '=' sign"));

    // Read 'zeros(n_rows,n_cols)' string
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Just read '" << tmpString << "'" << std::endl;
    unsigned int posInTmpString = 6;

    // Isolate 'n_rows' in a string
    char nRowsString[tmpString.size()-posInTmpString+1];
    unsigned int posInRowsString = 0;
    do {
      queso_require_less_msg(posInTmpString, tmpString.size(), "symbol ',' not found in first line of file");
      nRowsString[posInRowsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ',');
    nRowsString[posInRowsString] = '\0';

    // Isolate 'n_cols' in a string
    posInTmpString++; // Avoid reading ',' char
    char nColsString[tmpString.size()-posInTmpString+1];
    unsigned int posInColsString = 0;
    do {
      queso_require_less_msg(posInTmpString, tmpString.size(), "symbol ')' not found in first line of file");
      nColsString[posInColsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ')');
    nColsString[posInColsString] = '\0';

    // Convert 'n_rows' and 'n_cols' strings to numbers
    unsigned int numRowsInFile = (unsigned int) strtod(nRowsString,NULL);
    unsigned int numColsInFile = (unsigned int) strtod(nColsString,NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GslMatrix::subReadContents()"
                              << ": fullRank "        << m_env.fullRank()
                              << ", numRowsInFile = " << numRowsInFile
                              << ", numColsInFile = " << numColsInFile
                              << ", nRowsLocal = "    << nRowsLocal
                              << ", nCols = "         << nCols
                              << std::endl;
    }

    // Check if [num of rows in file] == [requested matrix row size]
    queso_require_equal_to_msg(numRowsInFile, nRowsLocal, "size of vec in file is not big enough");

    // Check if [num of cols in file] == [num cols in current matrix]
    queso_require_equal_to_msg(numColsInFile, nCols, "number of parameters of vec in file is different than number of parameters in this vec object");

    // Code common to any core in a communicator
    unsigned int maxCharsPerLine = 64*nCols; // Up to about 60 characters to represent each parameter value

    unsigned int lineId = 0;
    while (lineId < idOfMyFirstLine) {
      filePtrSet.ifsVar->ignore(maxCharsPerLine,'\n');
      lineId++;
    };

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GslMatrix::subReadContents()"
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
    queso_require_equal_to_msg(tmpString, std::string("="), std::string("in core 0, string should be the '=' sign"));

    // Take into account the ' [' portion
    std::streampos tmpPos = filePtrSet.ifsVar->tellg();
    filePtrSet.ifsVar->seekg(tmpPos+(std::streampos)2);

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GslMatrix::subReadContents()"
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
operator<<(std::ostream& os, const GslMatrix& obj)
{
  obj.print(os);

  return os;
}

GslMatrix operator*(double a, const GslMatrix& mat)
{
  GslMatrix answer(mat);
  answer *= a;
  return answer;
}

GslVector operator*(const GslMatrix& mat, const GslVector& vec)
{
  return mat.multiply(vec);
}

GslMatrix operator*(const GslMatrix& m1, const GslMatrix& m2)
{
  unsigned int m1Rows = m1.numRowsLocal();
  unsigned int m1Cols = m1.numCols();
  unsigned int m2Rows = m2.numRowsLocal();
  unsigned int m2Cols = m2.numCols();

  queso_require_equal_to_msg(m1Cols, m2Rows, "different sizes m1Cols and m2Rows");

  GslMatrix mat(m1.env(),m1.map(),m2Cols);

  //std::cout << "In GslMatrix(mat * mat): m1Cols = " << m1Cols << std::endl;

  unsigned int commonSize = m1Cols;
  for (unsigned int row1 = 0; row1 < m1Rows; ++row1) {
    for (unsigned int col2 = 0; col2 < m2Cols; ++col2) {
      double result = 0.;
      for (unsigned int k = 0; k < commonSize; ++k) {
        result += m1(row1,k)*m2(k,col2);
      }
      mat(row1,col2) = result;
    }
    //std::cout << "In GslMatrix(mat * mat): ended row " << row1 << std::endl;
  }

  return mat;
}

GslMatrix operator+(const GslMatrix& m1, const GslMatrix& m2)
{
  GslMatrix answer(m1);
  answer += m2;
  return answer;
}

GslMatrix operator-(const GslMatrix& m1, const GslMatrix& m2)
{
  GslMatrix answer(m1);
  answer -= m2;
  return answer;
}

GslMatrix matrixProduct(const GslVector& v1, const GslVector& v2)
{
  unsigned int nRows = v1.sizeLocal();
  unsigned int nCols = v2.sizeLocal();
  GslMatrix answer(v1.env(),v1.map(),nCols);

  for (unsigned int i = 0; i < nRows; ++i) {
    double value1 = v1[i];
    for (unsigned int j = 0; j < nCols; ++j) {
      answer(i,j) = value1*v2[j];
    }
  }

  return answer;
}

GslMatrix leftDiagScaling(const GslVector& vec, const GslMatrix& mat)
{
  unsigned int vSize = vec.sizeLocal();
  unsigned int mRows = mat.numRowsLocal();
  unsigned int mCols = mat.numCols();

  queso_require_equal_to_msg(vSize, mRows, "size of vector is different from the number of rows in matrix");

  queso_require_equal_to_msg(mCols, mRows, "routine currently works for square matrices only");

  GslMatrix answer(mat);
  for (unsigned int i = 0; i < mRows; ++i) {
    double vecValue = vec[i];
    for (unsigned int j = 0; j < mCols; ++j) {
      answer(i,j) *= vecValue;
    }
  }

  return answer;
}

GslMatrix rightDiagScaling(const GslMatrix& mat, const GslVector& vec)
{
  unsigned int vSize = vec.sizeLocal();
  unsigned int mRows = mat.numRowsLocal();
  unsigned int mCols = mat.numCols();

  queso_require_equal_to_msg(vSize, mCols, "size of vector is different from the number of cols in matrix");

  queso_require_equal_to_msg(mCols, mRows, "routine currently works for square matrices only");

  GslMatrix answer(mat);
  for (unsigned int j = 0; j < mCols; ++j) {
    double vecValue = vec[j];
    for (unsigned int i = 0; i < mRows; ++i) {
      answer(i,j) *= vecValue;
    }
  }

  return answer;
}

}  // End namespace QUESO
