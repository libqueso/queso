/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <uqGslMatrix.h>
#include <uqGslVector.h>
#include <uqDefines.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

uqGslMatrixClass::uqGslMatrixClass()
  :
  uqMatrixClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqGslMatrixClass::constructor(), default",
                      "should not be used by user");
}

uqGslMatrixClass::uqGslMatrixClass(
  const uqBaseEnvironmentClass& env,
  const Epetra_Map&             map,
  unsigned int                  numCols)
  :
  uqMatrixClass(env,map),
  m_mat        (gsl_matrix_calloc(map.NumGlobalElements(),numCols)),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.fullRank(),
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
                      m_env.fullRank(),
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
  m_mat        (gsl_matrix_calloc(v.sizeLocal(),v.sizeLocal())),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.fullRank(),
                      "uqGslMatrixClass::constructor(), eye",
                      "null matrix generated");

  for (unsigned int i = 0; i < m_mat->size1; ++i) {
    (*this)(i,i) = diagValue;
  }
}

uqGslMatrixClass::uqGslMatrixClass(const uqGslVectorClass& v)
  :
  uqMatrixClass(v.env(),v.map()),
  m_mat        (gsl_matrix_calloc(v.sizeLocal(),v.sizeLocal())),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.fullRank(),
                      "uqGslMatrixClass::constructor(), from vector",
                      "null matrix generated");

  unsigned int dim = v.sizeLocal();
  for (unsigned int i = 0; i < dim; ++i) {
    (*this)(i,i) = v[i];
  }
}

uqGslMatrixClass::uqGslMatrixClass(const uqGslMatrixClass& B)
  :
  uqMatrixClass(B.env(),B.map()),
  m_mat        (gsl_matrix_calloc(B.numRowsLocal(),B.numCols())),
  m_LU         (NULL),
  m_permutation(NULL)
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.fullRank(),
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
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  this->copy(obj);
  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator*=(double a)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  int iRC;
  iRC = gsl_matrix_scale(m_mat,a);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslMatrixClass::operator*=()",
                    "scaling failed");
  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator/=(double a)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  *this *= (1./a);

  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator+=(const uqGslMatrixClass& rhs)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  int iRC;
  iRC = gsl_matrix_add(m_mat,rhs.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslMatrixClass::operator+=()",
                    "failed");

  return *this;
}

uqGslMatrixClass&
uqGslMatrixClass::operator-=(const uqGslMatrixClass& rhs)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  int iRC;
  iRC = gsl_matrix_sub(m_mat,rhs.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslMatrixClass::operator-=()",
                    "failed");

  return *this;
}

double&
uqGslMatrixClass::operator()(unsigned int i, unsigned int j)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  UQ_FATAL_TEST_MACRO(i >= m_mat->size1,
                      m_env.fullRank(),
                      "uqGslMatrixClass::operator(i,j)",
                      "i is too large");
  UQ_FATAL_TEST_MACRO(j >= m_mat->size2,
                      m_env.fullRank(),
                      "uqGslMatrixClass::operator(i,j)",
                      "j is too large");
  return *gsl_matrix_ptr(m_mat,i,j);
}

const double&
uqGslMatrixClass::operator()(unsigned int i, unsigned int j) const
{
  UQ_FATAL_TEST_MACRO(i >= m_mat->size1,
                      m_env.fullRank(),
                      "uqGslMatrixClass::operator(i,j) const",
                      "i is too large");
  UQ_FATAL_TEST_MACRO(j >= m_mat->size2,
                      m_env.fullRank(),
                      "uqGslMatrixClass::operator(i,j) const",
                      "j is too large");
  return *gsl_matrix_const_ptr(m_mat,i,j);
}

void
uqGslMatrixClass::copy(const uqGslMatrixClass& src)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  int iRC;
  iRC = gsl_matrix_memcpy(this->m_mat, src.m_mat);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslMatrixClass::copy()",
                    "failed");

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

  unsigned int dim = this->numRowsLocal();
  double aux = 0.;
  for (unsigned int i = 0; i < dim; i++) {
    for (unsigned int j = 0; j < dim; j++) {
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

  unsigned int dim = this->numRowsLocal();
  double aux = 0.;
  for (unsigned int i = 0; i < dim; i++) {
    for (unsigned int j = 0; j < dim; j++) {
      aux = fabs((*this)(i,j));
      if (aux > value) {
        value = aux;
      }
    }
  }

  return value;
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
              m_env.fullRank(),
              "uqGslMatrixClass::chol()",
              "matrix is not positive definite",
              UQ_MATRIX_IS_NOT_POS_DEFINITE_RC);

  return iRC;
}

void
uqGslMatrixClass::zeroLower(bool includeDiagonal)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  if (this->numRowsLocal() != this->numCols()) return;

  unsigned int dim = this->numRowsLocal();
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
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  if (this->numRowsLocal() != this->numCols()) return;

  unsigned int dim = this->numRowsLocal();
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
  unsigned int nRows = this->numRowsLocal();
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
  UQ_FATAL_TEST_MACRO((this->numCols() != x.sizeLocal()),
                      m_env.fullRank(),
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
                      m_env.fullRank(),
                      "uqGslMatrixClass::multiply(), vector return void",
                      "matrix and x have incompatible sizes");

  UQ_FATAL_TEST_MACRO((y.sizeLocal() != x.sizeLocal()),
                      m_env.fullRank(),
                      "uqGslMatrixClass::multiply(), vector return void",
                      "y and x have incompatible sizes");

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
                      m_env.fullRank(),
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
                      m_env.fullRank(),
                      "uqGslMatrixClass::multiply(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.sizeLocal() != b.sizeLocal()),
                      m_env.fullRank(),
                      "uqGslMatrixClass::multiply(), return void",
                      "solution and rhs have incompatible sizes");

  int iRC;
  if (m_LU == NULL) {
    m_LU = gsl_matrix_calloc(numRowsLocal(),numCols());
    UQ_FATAL_TEST_MACRO((m_LU == NULL),
                        m_env.fullRank(),
                        "uqGslMatrixClass::invertMultiply()",
                        "gsl_matrix_calloc() failed");

    iRC = gsl_matrix_memcpy(m_LU, m_mat);
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.fullRank(),
                      "uqGslMatrixClass::invertMultiply()",
                      "gsl_matrix_memcpy() failed");

    m_permutation = gsl_permutation_calloc(numCols());
    UQ_FATAL_TEST_MACRO((m_permutation == NULL),
                        m_env.fullRank(),
                        "uqGslMatrixClass::invertMultiply()",
                        "gsl_permutation_calloc() failed");

    if (m_inDebugMode) {
      std::cout << "In uqGslMatrixClass::invertMultiply()"
                << ": before LU decomposition, m_LU = ";
      gsl_matrix_fprintf(stdout, m_LU, "%f");
      std::cout << std::endl;
    }

    int signum;
    iRC = gsl_linalg_LU_decomp(m_LU,m_permutation,&signum); 
    UQ_FATAL_RC_MACRO(iRC,
                      m_env.fullRank(),
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
                    m_env.fullRank(),
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

void              
uqGslMatrixClass::invertMultiply(const uqGslMatrixClass& B, uqGslMatrixClass& X) const
{
  
  // Sanity Checks
  UQ_FATAL_RC_MACRO( ((B.numRowsLocal() != X.numRowsLocal()) ||
		      (B.numCols()      != B.numCols()     )   ),
		      m_env.fullRank(),
		      "uqGslMatrixClass::invertMultiply()",
		      "Matrices B and X are incompatible");

  
  UQ_FATAL_RC_MACRO( ((this->numRowsLocal() != X.numRowsLocal()) ||
		      (this->numCols()      != X.numCols()     )   ),
		      m_env.fullRank(),
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
                      m_env.fullRank(),
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
                      m_env.fullRank(),
                      "uqGslMatrixClass::multiply(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.sizeLocal() != b.sizeLocal()),
                      m_env.fullRank(),
                      "uqGslMatrixClass::multiply(), return void",
                      "solution and rhs have incompatible sizes");

  int iRC;

  if( m_LU == NULL ) m_LU = gsl_matrix_calloc(numRowsLocal(),numCols());
  UQ_FATAL_TEST_MACRO((m_LU == NULL),
		      m_env.fullRank(),
		      "uqGslMatrixClass::invertMultiply()",
		      "gsl_matrix_calloc() failed");
  
  iRC = gsl_matrix_memcpy(m_LU, m_mat);
  UQ_FATAL_RC_MACRO(iRC,
		    m_env.fullRank(),
		    "uqGslMatrixClass::invertMultiply()",
		    "gsl_matrix_memcpy() failed");
  
  if( m_permutation == NULL ) m_permutation = gsl_permutation_calloc(numCols());
  UQ_FATAL_TEST_MACRO((m_permutation == NULL),
		      m_env.fullRank(),
		      "uqGslMatrixClass::invertMultiply()",
		      "gsl_permutation_calloc() failed");
  
  int signum;
  iRC = gsl_linalg_LU_decomp(m_LU,m_permutation,&signum); 
  UQ_FATAL_RC_MACRO(iRC,
		    m_env.fullRank(),
		    "uqGslMatrixClass::invertMultiply()",
		    "gsl_linalg_LU_decomp() failed");

  iRC = gsl_linalg_LU_solve(m_LU,m_permutation,b.data(),x.data()); 
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslMatrixClass::invertMultiply()",
                    "gsl_linalg_LU_solve() failed");

  return;
}

void
uqGslMatrixClass::eigen(uqGslVectorClass& eigenValues, uqGslMatrixClass* eigenVectors) const
{
  unsigned int n = eigenValues.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::eigen()",
                      "invalid input vector size");

  if (eigenVectors) {
    UQ_FATAL_TEST_MACRO((eigenValues.sizeLocal() != eigenVectors->numRowsLocal()),
                        UQ_UNAVAILABLE_RANK,
                        "uqGslVectorClass::eigen()",
                        "different input vector sizes");
  }

  if (eigenVectors == NULL) {
    gsl_eigen_symm_workspace* w = gsl_eigen_symm_alloc((const size_t) n);
    gsl_eigen_symm(m_mat,eigenValues.data(),w);
    gsl_eigen_symm_free(w);
  }
  else {
    gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc((const size_t) n);
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
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::largestEigen()",
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
                      UQ_UNAVAILABLE_RANK,
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
                      UQ_UNAVAILABLE_RANK,
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
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::smallestEigen()",
                      "Maximum num iterations exceeded");

  return;
}

void
uqGslMatrixClass::getColumn( const unsigned int column_num, uqGslVectorClass& column) const
{
  // Sanity checks
  UQ_FATAL_TEST_MACRO(( (column_num >= this->numCols()) ||
			(column_num < 0) ),
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::getColumn",
                      "Specified row number not within range");

  // Temporary working vector
  gsl_vector* gsl_column = gsl_vector_alloc( column.sizeLocal() );

  int error_code = gsl_matrix_get_col( gsl_column, m_mat, column_num );
  UQ_FATAL_RC_MACRO( (error_code != 0),
		     UQ_UNAVAILABLE_RANK,
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
uqGslMatrixClass::getRow( const unsigned int row_num, uqGslVectorClass& row) const
{
  // Sanity checks
  UQ_FATAL_TEST_MACRO(( (row_num >= this->numRowsLocal()) ||
			(row_num < 0) ),
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::getRow",
                      "Specified row number not within range");

  // Temporary working vector
  gsl_vector* gsl_row = gsl_vector_alloc( row.sizeLocal() );

  int error_code = gsl_matrix_get_row( gsl_row, m_mat, row_num );
  UQ_FATAL_RC_MACRO( (error_code != 0),
		     UQ_UNAVAILABLE_RANK,
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
uqGslMatrixClass::getRow( const unsigned int row_num ) const
{
  // We rely on the void getRow's to do sanity checks
  // So we don't do them here.

  uqGslVectorClass row(m_env, m_map);

  this->getRow( row_num, row );
  
  return row;
}

uqGslVectorClass
uqGslMatrixClass::getColumn( const unsigned int column_num ) const
{
  // We rely on the void getRow's to do sanity checks
  // So we don't do them here.

  uqGslVectorClass column(m_env, m_map);

  this->getColumn( column_num, column );
  
  return column;
}

void
uqGslMatrixClass::setRow( const unsigned int row_num, const uqGslVectorClass& row)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  // Sanity checks
  UQ_FATAL_TEST_MACRO(( (row_num >= this->numRowsLocal()) ||
			(row_num < 0) ),
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::setRow",
                      "Specified row number not within range");

  UQ_FATAL_TEST_MACRO((row.sizeLocal() != this->numCols()),
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::setRow",
                      "row vector not same size as this matrix");

  gsl_vector* gsl_row = row.data();

  int error_code = gsl_matrix_set_row( m_mat, row_num, gsl_row );
  UQ_FATAL_RC_MACRO( (error_code != 0),
		     UQ_UNAVAILABLE_RANK,
		     "uqGslMatrixClass::setRow()",
		     "gsl_matrix_set_row failed");

  return;
}

void
uqGslMatrixClass::setColumn( const unsigned int column_num, const uqGslVectorClass& column)
{
  if (m_LU) {
    gsl_matrix_free(m_LU);
    m_LU = NULL;
  }
  // Sanity checks
  UQ_FATAL_TEST_MACRO(( (column_num >= this->numCols()) ||
			(column_num < 0) ),
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::setColumn",
                      "Specified column number not within range");

  UQ_FATAL_TEST_MACRO((column.sizeLocal() != this->numRowsLocal()),
                      UQ_UNAVAILABLE_RANK,
                      "uqGslMatrixClass::setColumn",
                      "column vector not same size as this matrix");

  gsl_vector* gsl_column = column.data();

  int error_code = gsl_matrix_set_col( m_mat, column_num, gsl_column );
  UQ_FATAL_RC_MACRO( (error_code != 0),
		     UQ_UNAVAILABLE_RANK,
		     "uqGslMatrixClass::setColumn()",
		     "gsl_matrix_set_col failed");

  return;
}

void
uqGslMatrixClass::mpiSum( const MPI_Comm& comm, uqGslMatrixClass& M_global ) const
{
  // Sanity Checks
  UQ_FATAL_RC_MACRO( ( (this->numRowsLocal() != M_global.numRowsLocal()) ||
		       (this->numCols() != M_global.numCols())             ),
		       UQ_UNAVAILABLE_RANK,
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

  int ierr = MPI_Allreduce( (void*) &local[0], (void*) &global[0],
			    size, MPI_DOUBLE, MPI_SUM, comm );
  UQ_FATAL_TEST_MACRO(ierr != MPI_SUCCESS,
                      m_env.fullRank(),
                      "uqGslMatrixClass::mpiSum()",
                      "failed MPI_Allreduce()");

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

void
uqGslMatrixClass::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.fullRank(),
                      "uqGslMatrixClass::subWriteContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.fullRank(),
                      "uqGslMatrixClass::subWriteContents()",
                      "implemented just for sequential vectors for now");

  std::ofstream* ofsVar = NULL;
  m_env.openOutputFile(fileName,
                       UQ_FILE_EXTENSION_FOR_MATLAB_FORMAT,
                       allowedSubEnvIds,
                       false,
                       ofsVar);

  if (ofsVar) {
    unsigned int nRows = this->numRowsLocal();
    unsigned int nCols = this->numCols();
    *ofsVar << varNamePrefix << "_sub" << m_env.subIdString() << " = zeros(" << nRows
            << ","                                                           << nCols
            << ");"
            << std::endl;
    *ofsVar << varNamePrefix << "_sub" << m_env.subIdString() << " = [";

    for (unsigned int i = 0; i < nRows; ++i) {
      for (unsigned int j = 0; j < nCols; ++j) {
        *ofsVar << (*this)(i,j)
                << " ";
      }
      *ofsVar << "\n";
    }
    *ofsVar << "];\n";
    //ofsVar->close();
    delete ofsVar;
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
                      m1.env().fullRank(),
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

uqGslMatrixClass operator-(const uqGslMatrixClass& m1, const uqGslMatrixClass& m2)
{
  uqGslMatrixClass answer(m1);
  answer -= m2;
  return answer;
}

uqGslMatrixClass matrixProduct(const uqGslVectorClass& v1, const uqGslVectorClass& v2)
{
  unsigned int numRowsLocal = v1.sizeLocal();
  unsigned int numCols = v2.sizeLocal();
  uqGslMatrixClass answer(v1.env(),v1.map(),numCols);

  for (unsigned int i = 0; i < numRowsLocal; ++i) {
    double value1 = v1[i];
    for (unsigned int j = 0; j < numCols; ++j) {
      answer(i,j) = value1*v2[j];
    }
  }

  return answer;
}

uqGslMatrixClass diagScaling(const uqGslVectorClass& vec, const uqGslMatrixClass& mat)
{
  unsigned int vSize = vec.sizeLocal();
  unsigned int mRows = mat.numRowsLocal();
  unsigned int mCols = mat.numCols();

  UQ_FATAL_TEST_MACRO((vSize != mRows),
                      mat.env().fullRank(),
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
