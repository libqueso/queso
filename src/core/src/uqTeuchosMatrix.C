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

#include <uqDefines.h>

#ifdef QUESO_HAS_TRILINOS
#include <uqTeuchosMatrix.h>
#include <uqTeuchosVector.h>
#endif
#include <sys/time.h>
#include <cmath>

#ifdef QUESO_HAS_TRILINOS

namespace QUESO {

using std:: cout;
using std:: endl;
// ---------------------------------------------------
// default constructor ------------------------------- 
uqTeuchosMatrixClass::uqTeuchosMatrixClass()
  :
  uqMatrixClass()
{
  //this part is never called.
  UQ_FATAL_TEST_MACRO(true,		   
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::constructor(), default",
                      "should not be used by user");                   
}

// ---------------------------------------------------
// shaped constructor (square/rectangular)------------ 
uqTeuchosMatrixClass::uqTeuchosMatrixClass( // can be a rectangular matrix
  const uqBaseEnvironmentClass& env,
  const uqMapClass&             map,
  unsigned int                  nCols)
  :
  uqMatrixClass  (env,map),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  v_pivoting     (NULL), 
  m_signum       (0),
  m_isSingular   (false)
{
  m_mat.shape(map.NumGlobalElements(),nCols);
  m_LU.shape(0,0);   
}

// ---------------------------------------------------
// shaped constructor (square) ----------------------- 
uqTeuchosMatrixClass::uqTeuchosMatrixClass( // square matrix
  const uqBaseEnvironmentClass& env,
  const uqMapClass&             map,
  double                        diagValue)
  :
  uqMatrixClass  (env,map),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  v_pivoting     (NULL),  
  m_signum       (0),
  m_isSingular   (false)
{
  m_mat.shape    (map.NumGlobalElements(),map.NumGlobalElements());
  m_LU.shape(0,0); 

  for (unsigned int i = 0; i < (unsigned int) m_mat.numRows(); ++i) {
    m_mat(i,i) = diagValue;
  }
}
// ---------------------------------------------------
// shaped constructor (diagonal) --------------------- 
// Kemelli tested on 12/05/12
uqTeuchosMatrixClass::uqTeuchosMatrixClass( // square matrix
  const uqTeuchosVectorClass& v,
  double                      diagValue)
  :
  uqMatrixClass  (v.env(),v.map()),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  v_pivoting     (NULL), 
  m_signum       (0),
  m_isSingular   (false)
{
 m_mat.shape    (v.sizeLocal(),v.sizeLocal());
 m_LU.shape(0,0); 
 
 for (unsigned int i = 0; i < (unsigned int) m_mat.numRows(); ++i) {
    m_mat(i,i) = diagValue;
  }
}

// ---------------------------------------------------
// shaped constructor (diagonal, from vector) --------
// Kemelli, 12/05/12 - tested
uqTeuchosMatrixClass::uqTeuchosMatrixClass(const uqTeuchosVectorClass& v) // square matrix
  :
  uqMatrixClass  (v.env(),v.map()),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  v_pivoting     (NULL),
  m_signum       (0),
  m_isSingular   (false)
{
  m_mat.shape    (v.sizeLocal(),v.sizeLocal());
  m_LU.shape(0,0); 
  
  unsigned int dim =  v.sizeLocal();

  for (unsigned int i = 0; i < dim; ++i) {
    m_mat(i,i) = v[i];
  }
}
// ---------------------------------------------------
// copy constructor----------------------------------- 
// Kemelli 12/05/12 - tested
uqTeuchosMatrixClass::uqTeuchosMatrixClass(const uqTeuchosMatrixClass& B) // can be a rectangular matrix
  :
  uqMatrixClass  (B.env(),B.map()),
  m_inverse      (NULL),
  m_svdColMap    (NULL),
  m_svdUmat      (NULL),
  m_svdSvec      (NULL),
  m_svdVmat      (NULL),
  m_svdVTmat     (NULL),
  m_determinant  (-INFINITY),
  m_lnDeterminant(-INFINITY),
  v_pivoting     (NULL), 
  m_signum       (0),
  m_isSingular   (false)
{
  m_mat.shape    (B.numRowsLocal(),B.numCols());
  m_LU.shape(0,0); 
   
  this->uqMatrixClass::copy(B);
  this->copy(B);
}

// ---------------------------------------------------
// destructor ---------------------------------------- 
uqTeuchosMatrixClass::~uqTeuchosMatrixClass()
{
  this->resetLU();
}


// ---------------------------------------------------
// Set methodos --------------------------------------

uqTeuchosMatrixClass& uqTeuchosMatrixClass::operator=(const uqTeuchosMatrixClass& obj)
{
  this->resetLU();
  this->copy(obj);
  return *this;
}

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
uqTeuchosMatrixClass& uqTeuchosMatrixClass::operator*=(double a)
{
  this->resetLU();
  int iRC;
  iRC = m_mat.scale(a);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
                    "uqTeuchosMatrixClass::operator*=()",
                    "scaling failed");
  return *this;
}

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
uqTeuchosMatrixClass& uqTeuchosMatrixClass::operator/=(double a)
{
  this->resetLU();
  m_mat.scale(1./a);
  return *this;
}

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
uqTeuchosMatrixClass& uqTeuchosMatrixClass::operator+=(const uqTeuchosMatrixClass& rhs)
{
  this->resetLU();

  unsigned int i,j, nrows=rhs.numRowsLocal(), ncols=rhs.numCols();
  
  for(i=0; i< nrows ; i++)
    for (j = 0; j < ncols; j++) 
      (*this)(i,j) += rhs(i,j);
    
  return *this;
}
// ---------------------------------------------------
// Kemelli 12/05/12 - tested
uqTeuchosMatrixClass&
uqTeuchosMatrixClass::operator-=(const uqTeuchosMatrixClass& rhs)
{
  this->resetLU();

  unsigned int i,j, nrows=rhs.numRowsLocal(), ncols=rhs.numCols();
  
  for(i=0; i< nrows ; i++)
    for (j = 0; j < ncols; j++) 
      (*this)(i,j) -= rhs(i,j);
    
  return *this;
}

// ---------------------------------------------------
// Accessor methods ----------------------------------

// Kemelli 12/05/12 - tested
double& uqTeuchosMatrixClass::operator()(unsigned int i, unsigned int j)
{
  this->resetLU();
  if ((i >= (unsigned int) m_mat.numRows()) ||
      (j >= (unsigned int) m_mat.numCols())) {
    std::cerr << "In uqTeuchosMatrixClass::operator(i,j) (non-const)"
              << ": i = " << i
              << ", j = " << j
              << ", m_mat.numRows() = " << (unsigned int) m_mat.numRows()
              << ", m_mat.numCols() = " << (unsigned int) m_mat.numCols()
              << std::endl;
    UQ_FATAL_TEST_MACRO(i >= (unsigned int) m_mat.numRows(),
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::operator(i,j)",
                        "i is too large");
    UQ_FATAL_TEST_MACRO(j >= (unsigned int) m_mat.numCols(),
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::operator(i,j)",
                        "j is too large");
  }
  return m_mat(i,j);
}

// ---------------------------------------------------
const double& uqTeuchosMatrixClass::operator()(unsigned int i, unsigned int j) const
{
  UQ_FATAL_TEST_MACRO(i >= (unsigned int) m_mat.numRows(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::operator(i,j) const",
                      "i is too large");
  UQ_FATAL_TEST_MACRO(j >= (unsigned int) m_mat.numCols(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::operator(i,j) const",
                      "j is too large");
  return m_mat(i,j);
}

// ---------------------------------------------------
// Attribute methods ---------------------------------

// Kemelli 12/05/12 - tested 
unsigned int
uqTeuchosMatrixClass::numRowsLocal() const
{
  return (unsigned int) m_mat.numRows();
}

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
unsigned int 
uqTeuchosMatrixClass::numRowsGlobal() const
{
  return (unsigned int) m_mat.numRows();
}

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
unsigned int 
uqTeuchosMatrixClass::numCols() const
{
  return (unsigned int) m_mat.numCols();
};

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
double*
uqTeuchosMatrixClass::values()
{
  return  m_mat.values();
};

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
int
uqTeuchosMatrixClass::stride()
{
  return  m_mat.stride();
};



// ---------------------------------------------------
double
uqTeuchosMatrixClass::max() const
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


// ---------------------------------------------------
unsigned int
uqTeuchosMatrixClass::rank(double absoluteZeroThreshold, double relativeZeroThreshold) const
{
  int iRC = 0;
  iRC = internalSvd();
  if (iRC) {}; // just to remove compiler warning

  uqTeuchosVectorClass relativeVec(*m_svdSvec);
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
    *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::rank()"
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

// ---------------------------------------------------
// checked 1/07/13
uqTeuchosMatrixClass
uqTeuchosMatrixClass::transpose() const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((nRows != nCols),
                       m_env.worldRank(),
                       "uqTeuchosMatrixClass::transpose()",
                       "routine works only for square matrices");

  uqTeuchosMatrixClass mat(m_env,m_map,nCols);
  for (unsigned int row = 0; row < nRows; ++row) {
    for (unsigned int col = 0; col < nCols; ++col) {
      mat(row,col) = (*this)(col,row);
    }
  }

  return mat;
}

// ---------------------------------------------------
// tested 1/31/13
uqTeuchosMatrixClass
uqTeuchosMatrixClass::inverse() const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((nRows != nCols),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::inverse()",
                      "matrix is not square");

  if (m_inverse == NULL) {
    m_inverse = new uqTeuchosMatrixClass(m_env,m_map,nCols);
    uqTeuchosVectorClass unitVector(m_env,m_map);
    unitVector.cwSet(0.);
    uqTeuchosVectorClass multVector(m_env,m_map);
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
    *m_env.subDisplayFile() << "CHECKING In uqTeuchosMatrixClass::inverse()"
                            << ": M.lnDet = "      << this->lnDeterminant()
                            << ", M^{-1}.lnDet = " << m_inverse->lnDeterminant()
                            << std::endl;
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::inverse():"
                            << "\n M = "        << *this
                            << "\n M^{-1} = "   << *m_inverse
                            << "\n M*M^{-1} = " << (*this)*(*m_inverse)
                            << "\n M^{-1}*M = " << (*m_inverse)*(*this)
                            << std::endl;
  }

/* cout << "Inverse ---------------------\n";
   for (int i=0; i<m_inverse->numRowsLocal(); i++) {
     for (int j=0; j<m_inverse->numCols(); j++)
       cout<< m_inverse->m_mat(i,j) << " ";
     cout << endl;
   }
   cout << endl;
*/
  return *m_inverse;
}

// ----------------------------------------------
//checked 12/10/12
double
uqTeuchosMatrixClass::determinant() const
{
  if (m_determinant == -INFINITY) 
  {    
    if(m_LU.numRows() ==0 && m_LU.numCols() ==0)  //dummy
    {
      uqTeuchosVectorClass tmpB(m_env,m_map);
      uqTeuchosVectorClass tmpX(m_env,m_map);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::lnDeterminant()"
                                << ": before 'this->invertMultiply()'"
                                << std::endl;
      }
      this->invertMultiply(tmpB,tmpX);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::lnDeterminant()"
                                << ": after 'this->invertMultiply()'"
                                << std::endl;
      }
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::lnDeterminant()"
                              << ": before computing det"
                              << std::endl;
    }
    
    double det   = 1.0;
    double lnDet = 0.0;
    for (int i=0;i<m_LU.numCols();i++) {
      det   *= m_LU(i,i);
      lnDet += std::log(m_LU(i,i));
    }
  
    m_determinant   = det;
    m_lnDeterminant = lnDet;
        
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::lnDeterminant()"
                              << ": after computing det"
                              << std::endl;
    }
  }

  return m_determinant;
}

// ----------------------------------------------
//checked 12/10/12
double
uqTeuchosMatrixClass::lnDeterminant() const
{
  if (m_lnDeterminant == -INFINITY) 
  {
    if(m_LU.numRows() ==0 && m_LU.numCols() ==0)  //dummy
    {
      uqTeuchosVectorClass tmpB(m_env,m_map);
      uqTeuchosVectorClass tmpX(m_env,m_map);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::lnDeterminant()"
                                << ": before 'this->invertMultiply()'"
                                << std::endl;
      }
      this->invertMultiply(tmpB,tmpX);
      if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
        *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::lnDeterminant()"
                                << ": after 'this->invertMultiply()'"
                                << std::endl;
      }
    }
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::lnDeterminant()"
                              << ": before computing lnDet"
                              << std::endl;
    }
    
    double det   = 1.0;
    double lnDet = 0.0;
    for (int i=0;i<m_LU.numCols();i++) {
      det   *= m_LU(i,i);
      lnDet += std::log(m_LU(i,i));
    }

    m_determinant   = det;
    m_lnDeterminant = lnDet;
    
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::lnDeterminant()"
                              << ": after computing lnDet"
                              << std::endl;
    }
  }

  return m_lnDeterminant;
}

// ---------------------------------------------------
// Norm methods --------------------------------------

double
uqTeuchosMatrixClass::normFrob() const
{
  return m_mat.normFrobenius ();
}

// ---------------------------------------------------
double
uqTeuchosMatrixClass::normMax() const
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

// Mathematical methods ------------------------------
// ---------------------------------------------------

// Kemelli: added and tested 12/10/12
int uqTeuchosMatrixClass::chol()
{
  int return_success =0 ;
/*  If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A (lets call it L), and the strictly upper
 *          triangular part of A is not referenced. Thus, the upper triangular part 
 *          of A +++must be manually+++ overwritten with L^T. */
  
  Teuchos::LAPACK<int, double> lapack;
  int info;//= 0:  successful exit
           //< 0:  if INFO = -i, the i-th argument had an illegal value
           //> 0:  if INFO = i, the leading minor of order i is not
           //      positive definite, and the factorization could not be
           //      completed.
  char UPLO = 'U'; 
          //= 'U':  Upper triangle of A is stored;
          //= 'L':  Lower triangle of A is stored.
  
  lapack.POTRF (UPLO, m_mat.numRows(), m_mat.values(), m_mat.stride(), &info);
  
  // Overwriting the upper triangular part of the input matrix A with  L^T 
  //(the diagonal terms are identical for both L and L^T)

  for (int i=0;i<m_mat.numRows();i++){
    for (int j=i+1;j<m_mat.numCols();j++)
       m_mat(i,j) = m_mat(j,i) ;      
  }
  
  if (info != 0) {
    std::cerr << "In uqTeuchosMtrixClass::chol()"
              << ": INFO = " << info
              << ",\n INFO < 0:  if INFO = -i, the i-th argument had an illegal value."
              << "\n INFO > 0:  if INFO = i, the leading minor of order i is not "
              << " positive definite, and the factorization could not be completed."
              << std::endl;
    return_success =1 ;
  }
  
  UQ_RC_MACRO(info, // Yes, *not* a fatal check on RC
              m_env.worldRank(),
              "uqTeuchosMatrixClass::chol()",
              "matrix is not positive definite",
              UQ_MATRIX_IS_NOT_POS_DEFINITE_RC);
     
  return return_success;  
};

// ---------------------------------------------------
int
uqTeuchosMatrixClass::svd(uqTeuchosMatrixClass& matU, uqTeuchosVectorClass& vecS, uqTeuchosMatrixClass& matVt) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((matU.numRowsLocal() != nRows) || (matU.numCols() != nCols),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::svd()",
                      "invalid matU");

  UQ_FATAL_TEST_MACRO((vecS.sizeLocal() != nCols), //std::min(nRows,nRows)),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::svd()",
                      "invalid vecS");

  UQ_FATAL_TEST_MACRO((matVt.numRowsLocal() != nCols) || (matVt.numCols() != nCols),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::svd()",
                      "invalid matVt");

  int iRC = internalSvd();

  matU  = *m_svdUmat;
  vecS  = *m_svdSvec;
  matVt = *m_svdVTmat;

  return iRC;
}

//---------------------------------------------------------------
//tested 2/1/13
const uqTeuchosMatrixClass& uqTeuchosMatrixClass::svdMatU() const
{
  int iRC = 0;
  iRC = internalSvd();
  if (iRC) {}; // just to remove compiler warning
  return *m_svdUmat;
}

//---------------------------------------------------------------
//tested 2/1/13
const uqTeuchosMatrixClass& uqTeuchosMatrixClass::svdMatV() const
{
  int iRC = 0;
  iRC = internalSvd();
  if (iRC) {}; // just to remove compiler warning
  return *m_svdVmat;
}

//---------------------------------------------------------------
// checked 2/27/13
/* An orthogonal matrix M has a norm-preserving property, i.e. 
 * for any vector v, \f[||Mv|| = ||v|| \f] (1). Then: 
 * \f[ min(||Ax − b||^2) = min(||Ax − b||) = min(||UDVT x − b||) = (1) min(||DV x − U b||) \f].
 * Substituting \f[ y = VT x \f] and \f[ b' = UT b \f] gives us \f[ Dy = b' \f] 
 * with D a diagonal matrix. Or,  \f[ y = inv(D)*UT*b \f] and we only have to 
 * solve the linear system: \f[ VT x = y \f].
 */
int
uqTeuchosMatrixClass::svdSolve(const uqTeuchosVectorClass& rhsVec, uqTeuchosVectorClass& solVec) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();
  unsigned int i;
  
  UQ_FATAL_TEST_MACRO((rhsVec.sizeLocal() != nRows),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::svdSolve()",
                      "invalid rhsVec");

  UQ_FATAL_TEST_MACRO((solVec.sizeLocal() != nCols),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::svdSolve()",
                      "invalid solVec");

  int iRC = internalSvd();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::svdSolve():"
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

 if (iRC == 0)
 {
   uqTeuchosMatrixClass  invD = uqTeuchosMatrixClass(solVec); 
   uqTeuchosMatrixClass  auxMatrix = uqTeuchosMatrixClass(solVec); 
    
   for (i=0; i<nRows; i++){
     invD(i,i) = 1./(m_svdSvec->values()[i]); 
   }
 
  // GESV: For a system Ax=b, on entry, b contains the right-hand side b; 
  // on exit it contains the solution x.
  // Thus, intead of doing y = inv(D)*UT*rhsVec = auxMatrix*rhsVec, with 
  // auxMatrix = inv(D)*UT; and then assigning solVec=y (requirement for using GESV)
  // lets do solVec= auxMatrix*rhsVec, and save one step in the calculation 
 
    auxMatrix = invD * svdMatU().transpose();
    solVec = auxMatrix*rhsVec;

  // solve the linear system VT * solVec = y
  // GESV changes the values of the matrix, so lets make a copy of it and use it.
  
    uqTeuchosMatrixClass* aux_m_svdVTmat  = new uqTeuchosMatrixClass(*m_svdVTmat);
   
    int ipiv[0], info;
    Teuchos::LAPACK<int, double> lapack;
    lapack.GESV(nCols, 1,  aux_m_svdVTmat->values(),  aux_m_svdVTmat->stride(), ipiv, &solVec[0], solVec.sizeLocal(), &info );  
  
  /* GESV output INFO: = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
   *                has been completed, but the factor U is exactly
   *                singular, so the solution could not be computed.   */
   
    iRC = info;
    delete aux_m_svdVTmat;    
  }
  return iRC; 
}

// ---------------------------------------------------
//checked 2/27/13
int
uqTeuchosMatrixClass::svdSolve(const uqTeuchosMatrixClass& rhsMat, uqTeuchosMatrixClass& solMat) const
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((rhsMat.numRowsLocal() != nRows),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::svdSolve()",
                      "invalid rhsMat");

  UQ_FATAL_TEST_MACRO((solMat.numRowsLocal() != nCols),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::svdSolve()",
                      "invalid solMat");

  UQ_FATAL_TEST_MACRO((rhsMat.numCols() != solMat.numCols()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::svdSolve()",
                      "rhsMat and solMat are not compatible");

  uqTeuchosVectorClass rhsVec(m_env,rhsMat.map());
  uqTeuchosVectorClass solVec(m_env,solMat.map());
  int iRC = 0;
  for (unsigned int j = 0; j < rhsMat.numCols(); ++j) {
    rhsVec = rhsMat.getColumn(j);
    iRC = this->svdSolve(rhsVec, solVec);
    if (iRC) break;
    solMat.setColumn(j,solVec);
  }

  return iRC;
}

// ---------------------------------------------------
// implemented/checked 1/8/13
// multiply this matrix by vector x and store in vector y, which is returned
// eg: uqTeuchosVectorClass v1(paramSpace.zeroVector() ); 
//     v1=covMatrix.multiply(meanVector);
uqTeuchosVectorClass
uqTeuchosMatrixClass::multiply(const uqTeuchosVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != x.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::multiply(), return vector",
                      "matrix and vector have incompatible sizes");

  uqTeuchosVectorClass y(m_env,m_map);
  this->multiply(x,y);

  return y;
}

// ---------------------------------------------------
//Kemelli checked 12/06/12
uqTeuchosVectorClass
uqTeuchosMatrixClass::invertMultiply(const uqTeuchosVectorClass& b) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply(), return vector",
                      "matrix and rhs have incompatible sizes");
  uqTeuchosVectorClass x(m_env,m_map);
  
  this->invertMultiply(b,x);

  return x;
}

// ---------------------------------------------------
//Kemelli checked 12/06/12
void
uqTeuchosMatrixClass::invertMultiply(const uqTeuchosVectorClass& b, uqTeuchosVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.sizeLocal() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply(), return void",
                      "solution and rhs have incompatible sizes");
 
  if (m_LU.numCols() == 0 && m_LU.numRows() == 0)  
  {						        
  UQ_FATAL_TEST_MACRO((v_pivoting != NULL),
                         m_env.worldRank(),
                         "uqTeuchosMatrixClass::invertMultiply()",
                         "v_pivoting should be NULL");
    
  //allocate m_LU and v_pivoting
  m_LU = m_mat;
  v_pivoting =(int *) malloc(sizeof(int)*m_LU.numCols() );
  
  UQ_FATAL_TEST_MACRO((m_LU.numCols() == 0 && m_LU.numRows() == 0),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply()",
                      "malloc() failed");

  UQ_FATAL_TEST_MACRO((v_pivoting == NULL),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply()",
                      "malloc() failed");

  if (m_inDebugMode) {
    std::cout << "In uqTeuchosMatrixClass::invertMultiply()"
	      << ": before LU decomposition, m_LU = ";
	      for (int i=0;i<3;i++){
		for (int j=0;j<3;j++)
		  cout << m_LU(i,j) <<"\t" ;
		cout << endl;
	      }
    }
      
  // Perform an LU factorization of matrix m_LU. Checked 12/06/12
  Teuchos::LAPACK<int, double> lapack;
  int info;
   
  lapack.GETRF( m_LU.numRows(), m_LU.numCols(), m_LU.values(), m_LU.stride(), v_pivoting, &info ); 
  
  if (info != 0) {
    std::cerr   << "In uqTeuchosMatrixClass::invertMultiply()"
		<< ", after lapack.GETRF"
                << ": INFO = " << info
                << ",\nINFO < 0:  if INFO = -i, the i-th argument had an illegal value.\n"
                << "INFO > 0:  if INFO = i, U(i,i) is exactly zero. The factorization \n" 
                << "has been completed, but the factor U is exactly singular, and division \n"
                << "by zero will occur if it is used to solve a system of equations."
                << std::endl;
  } 
  UQ_FATAL_RC_MACRO(info,
		    m_env.worldRank(),
		    "uqTeuchosMatrixClass::invertMultiplyForceLU()",
		    "GETRF() failed");
    
  if (info >  0) 
     m_isSingular = true;
     
  if (m_inDebugMode) {
    std::cout << "In uqTeuchosMatrixClass::invertMultiply()"
              << ": after LU decomposition, m_LU = ";
	      for (int i=0;i<3;i++){
		for (int j=0;j<3;j++)
		  cout << m_LU(i,j) <<"\t" ;
		cout << endl;
	      }      
      std::cout << std::endl;
    }       
  }
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::invertMultiply()"
                            << ": before 'lapack.GETRS()'"
                            << std::endl;
  }
  
  // Solve the linear system.
  Teuchos::LAPACK<int, double> lapack;
  int NRHS = 1; // NRHS: number of right hand sides, i.e., the number 
				// of columns of the matrix B. In this case, vector b.
  char TRANS = 'N';  // 'N':  A * x= B  (No transpose). Specifies the 
				     // form of the system of equations.
  int info02;
  
  //GETRS expects the matrix to be already factored in LU and uses the 
  //same ipiv vector, which are the pivot indices of the LU factorization

  x=b;                
  lapack.GETRS(TRANS, m_LU.numRows(), NRHS, m_LU.values(), m_LU.stride(), v_pivoting, &x[0],x.sizeLocal(), &info02 );
 
  if (info02 != 0) {
      std::cerr << "In uqTeuchosMatrixClass::invertMultiply()"
                << ", after lapack.GETRS - solve LU system"
                << ": INFO = " << info02
                << ",\nINFO < 0:  if INFO = -i, the i-th argument had an illegal value.\n"
                << std::endl;
  } 
  UQ_FATAL_RC_MACRO(info02,
		    m_env.worldRank(),
		    "uqTeuchosMatrixClass::invertMultiplyForceLU()",
		    "GETRS() failed"); 
 if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
    *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::invertMultiply()"
			    << ", after lapack.GETRS() - solve LU system."
                            << std::endl;
 }
 if (m_inDebugMode) {
    uqTeuchosVectorClass tmpVec(b - (*this)*x);
    std::cout << "In uqTeuchosMatrixClass::invertMultiply()"
              << ": ||b - Ax||_2 = "         << tmpVec.norm2()
              << ": ||b - Ax||_2/||b||_2 = " << tmpVec.norm2()/b.norm2()
              << std::endl;
 }
 return;
}

// ----------------------------------------------
uqTeuchosMatrixClass
uqTeuchosMatrixClass::invertMultiply(const uqTeuchosMatrixClass& B) const
{
  uqTeuchosMatrixClass X(m_env,m_map,B.numCols());
  this->invertMultiply(B,X);
  return X;
}

// ----------------------------------------------
//checked 12/10/12
void              
uqTeuchosMatrixClass::invertMultiply(const uqTeuchosMatrixClass& B, uqTeuchosMatrixClass& X) const
{
  // Sanity Checks
  UQ_FATAL_RC_MACRO(((B.numRowsLocal() != X.numRowsLocal()) || (B.numCols() != X.numCols())),
             m_env.worldRank(),
		    "uqTeuchosMatrixClass::invertMultiply()",
		    "Matrices B and X are incompatible");
  
  UQ_FATAL_RC_MACRO((this->numRowsLocal() != X.numRowsLocal()),
             m_env.worldRank(),
		    "uqTeuchosMatrixClass::invertMultiply()",
		    "This and X matrices are incompatible");

  // Some local variables used within the loop.
  uqTeuchosVectorClass b(m_env, m_map);
  uqTeuchosVectorClass x(m_env, m_map);
  
  for( unsigned int j = 0; j < B.numCols(); ++j ) {
    b = B.getColumn(j);
    //invertMultiply will only do the LU factorization once and store it. 
    x = this->invertMultiply(b);
    X.setColumn(j,x);
  }
  return;
}

//-----------------------------------------------
uqTeuchosVectorClass
uqTeuchosMatrixClass::invertMultiplyForceLU(const uqTeuchosVectorClass& b) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply(), return vector",
                      "matrix and rhs have incompatible sizes");

  uqTeuchosVectorClass x(m_env,m_map);
  this->invertMultiplyForceLU(b,x);

  return x;
}

// ---------------------------------------------------
// implemented and checked on 1/9/13
void
uqTeuchosMatrixClass::invertMultiplyForceLU(const uqTeuchosVectorClass& b, uqTeuchosVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.sizeLocal() != b.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply(), return void",
                      "solution and rhs have incompatible sizes");
 
  if (m_LU.numCols() == 0 && m_LU.numRows() == 0) {						      
    UQ_FATAL_TEST_MACRO((v_pivoting != NULL),
                         m_env.worldRank(),
                         "uqTeuchosMatrixClass::invertMultiply()",
                         "v_pivoting should be NULL");   
  }
  
  //allocate m_LU , yes outside the if above
  m_LU = m_mat;
  
  UQ_FATAL_TEST_MACRO((m_LU.numCols() == 0 && m_LU.numRows() == 0),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply()",
                      "Teuchos atttribuition m_LU = m_mat failed");

  //allocate v_pivoting
  if ( v_pivoting == NULL )
    v_pivoting =(int *) malloc(sizeof(int)*m_LU.numCols() );

  UQ_FATAL_TEST_MACRO((v_pivoting == NULL),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::invertMultiply()",
                      "malloc() for v_pivoting failed");
      
  // Perform an LU factorization of matrix m_LU. Checked 12/06/12
  Teuchos::LAPACK<int, double> lapack;
  int info;
   
  lapack.GETRF( m_LU.numRows(), m_LU.numCols(), m_LU.values(), m_LU.stride(), v_pivoting, &info ); 
  
  if (info != 0) {
    std::cerr   << "In uqTeuchosMatrixClass::invertMultiply()"
		<< ", after lapack.GETRF"
                << ": INFO = " << info
                << ",\nINFO < 0:  if INFO = -i, the i-th argument had an illegal value.\n"
                << "INFO > 0:  if INFO = i, U(i,i) is exactly zero. The factorization \n" 
                << "has been completed, but the factor U is exactly singular, and division \n"
                << "by zero will occur if it is used to solve a system of equations."
                << std::endl;
  } 
  
  UQ_FATAL_RC_MACRO(info,
		    m_env.worldRank(),
		    "uqTeuchosMatrixClass::invertMultiplyForceLU()",
		    "GETRF() failed");
 
  if (info >  0) 
     m_isSingular = true;
  
  // Solve the linear system.
  int NRHS = 1; // NRHS: number of right hand sides, i.e., the number 
				// of columns of the matrix B. In this case, vector b.
  char TRANS = 'N';  // 'N':  A * x= B  (No transpose). Specifies the 
				     // form of the system of equations.
  int info02;
  
  //GETRS expects the matrix to be already factored in LU and uses the 
  //same ipiv vector, which are the pivot indices of the LU factorization

  x=b;                
  lapack.GETRS(TRANS, m_LU.numRows(), NRHS, m_LU.values(), m_LU.stride(), v_pivoting, &x[0],x.sizeLocal(), &info02 );
 
  if (info02 != 0) {
      std::cerr << "In uqTeuchosMatrixClass::invertMultiply()"
                << ", after lapack.GETRS - solve LU system"
                << ": INFO = " << info02
                << ",\nINFO < 0:  if INFO = -i, the i-th argument had an illegal value.\n"
                << std::endl;
    }  
    
    UQ_FATAL_RC_MACRO(info02,
		    m_env.worldRank(),
		    "uqTeuchosMatrixClass::invertMultiplyForceLU()",
		    "GETRS() failed"); 
   
 return;
}
// ---------------------------------------------------
// Implemented and checked on 1/7/2013

/*  DSYEV computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.
*
*  Arguments
*  JOBZ    = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*  UPLO    = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*  N       The order of the matrix A.  N >= 0.
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N).On entry, the symmetric 
*  matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of A contains the
*  upper triangular part of the matrix A.  If UPLO = 'L', the leading N-by-N lower triangular 
*  part of A contains the lower triangular part of the matrix A. On exit, if JOBZ = 'V', then
*  if INFO = 0, A contains the orthonormal eigenvectors of the matrix A. If JOBZ = 'N', then 
*  on exit the lower triangle (if UPLO='L') or the upper triangle (if UPLO='U') of A, 
*  including the  diagonal, is destroyed.
*  LDA     (input) INTEGER - the leading dimension of the array A.  LDA >= max(1,N).
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*  LWORK   (input) INTEGER - the length of the array WORK.  LWORK >= max(1,3*N-1).
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal elements of 
*                an intermediate tridiagonal form did not converge to zero.*/
void
uqTeuchosMatrixClass::eigen(uqTeuchosVectorClass& eigenValues, uqTeuchosMatrixClass* eigenVectors) const
{
  unsigned int n = eigenValues.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::eigen()",
                      "invalid input vector size");

  if (eigenVectors) {
    UQ_FATAL_TEST_MACRO((eigenValues.sizeLocal() != eigenVectors->numRowsLocal()),
                        env().fullRank(),
                        "uqTeuchosVectorClass::eigen()",
                        "different input vector sizes");
  }
  
  // At the end of execution, Lapack funcion destroys the input matrix
  // So, lets not use m_mat, but instead, a copy of it
  Teuchos::SerialDenseMatrix<int,double> copy_m_mat;
  copy_m_mat = m_mat;

  Teuchos::LAPACK<int, double> lapack;
  int lwork = 3*n -1; 
  int info;
  char UPLO = 'L'; 
  double *W, *WORK;

  W = new double[n];
  WORK = new double[lwork];
      
  if (eigenVectors == NULL) {
    char JOBZ = 'N';//eigenvalues only        
     
    lapack.SYEV(JOBZ, UPLO, copy_m_mat.numRows(), copy_m_mat.values(), copy_m_mat.stride(), &W[0], &WORK[0], lwork, &info);
  
    for (unsigned int i=0; i< n; i++) 
      eigenValues[i] = W[i];
  
    UQ_FATAL_TEST_MACRO((info != 0),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::eigen()",
                      "invalid input vector size");
  }
  else {
  	char JOBZ = 'V';//eigenvalues and eigenvectors          
    
    lapack.SYEV(JOBZ, UPLO, copy_m_mat.numRows(), copy_m_mat.values(), copy_m_mat.stride(), &W[0], &WORK[0], lwork, &info);
  
  	for (unsigned int i=0; i< n; i++)
    	eigenValues[i] = W[i];
  
  	eigenVectors->m_mat = copy_m_mat;
  } 

  if (info != 0) {
    std::cerr << "In uqTeuchosMtrixClass::eigen()"
              << ": INFO = " << info
              << ",\n INFO < 0:  if INFO = -i, the i-th argument had an illegal value."
              << "\n INFO > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal "
		      << " elements of an intermediate tridiagonal form did not converge to zero."
              << std::endl;     
  }

  return;
}
// ---------------------------------------------------
void
uqTeuchosMatrixClass::largestEigen(double& eigenValue, uqTeuchosVectorClass& eigenVector) const
{ 
  unsigned int n = eigenVector.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::largestEigen()",
                      "invalid input vector size");
  Teuchos::LAPACK<int, double> lapack;
  
  //SYEV destroys the input matrix, so lets operate in a copy  
  Teuchos::SerialDenseMatrix<int,double> copy_m_mat;
  copy_m_mat = m_mat;

  int lwork = 3*n -1; 
  int info;
  char UPLO = 'L'; 
  char JOBZ = 'V'; //eigenvalues and eigenvectors          
  double *W, *WORK;
  
  W = new double[n];
  WORK = new double[lwork];
      
  lapack.SYEV(JOBZ, UPLO, copy_m_mat.numRows(), copy_m_mat.values(), copy_m_mat.stride(), &W[0], &WORK[0], lwork, &info);
   
  if (info != 0) {
    std::cerr << "In uqTeuchosMtrixClass::largestEigen()"
              << ": INFO = " << info
              << ",\n INFO < 0:  if INFO = -i, the i-th argument had an illegal value."
              << "\n INFO > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal "
		      << " elements of an intermediate tridiagonal form did not converge to zero."
              << std::endl;    
  }
  // If INFO = 0, W contains the eigenvalues in ascending order. 
  // Thus the largest eigenvalue is in W[n-1].
  eigenValue = W[n-1];
  
  // Eigenvector associated to the largest eigenvalue.
  // Stored in the n-th column of matrix copy_m_mat.
  for (int i=0; i< copy_m_mat.numRows(); i++)
    eigenVector[i] = copy_m_mat(i,n-1);
   
  return;
}

// ---------------------------------------------------
void
uqTeuchosMatrixClass::smallestEigen(double& eigenValue, uqTeuchosVectorClass& eigenVector) const
{ 
  unsigned int n = eigenVector.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::smallestEigen()",
                      "invalid input vector size");
  Teuchos::LAPACK<int, double> lapack;
  
  //SYEV destroys the input matrix, so lets operate in a copy  
  Teuchos::SerialDenseMatrix<int,double> copy_m_mat;
  copy_m_mat = m_mat;

  int lwork = 3*n -1; 
  int info;
  char UPLO = 'L'; 
  char JOBZ = 'V';//eigenvalues and eigenvectors          
  double *W, *WORK;
  
  W = new double[n];
  WORK = new double[lwork];
     
  lapack.SYEV(JOBZ, UPLO, copy_m_mat.numRows(), copy_m_mat.values(), copy_m_mat.stride(), &W[0], &WORK[0], lwork, &info);
  
  if (info != 0) {
    std::cerr << "In uqTeuchosMtrixClass::smallestEigen()"
              << ": INFO = " << info
              << ",\n INFO < 0:  if INFO = -i, the i-th argument had an illegal value."
              << "\n INFO > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal "
		      << " elements of an intermediate tridiagonal form did not converge to zero."
              << std::endl;    
  }
  
  // If INFO = 0, W contains the eigenvalues in ascending order. 
  // Thus the smallest eigenvalue is in W[0].
  eigenValue = W[0];
  
  // Eigenvector associated to the smallest eigenvalue.
  // Stored in the n-th column of matrix copy_m_mat.
  for (int i=0; i< copy_m_mat.numRows(); i++)
    eigenVector[i] = copy_m_mat(i,0);
   
  return;
}

// Get/Set methodos ----------------------------------
// ---------------------------------------------------

// ---------------------------------------------------
void
uqTeuchosMatrixClass::cwSet(double value)
{
  m_mat.putScalar(value);
  return;
}

// ---------------------------------------------------
// tested on 1/31/13: copies matrix mat to this matrix, starting at row rowId and column colId
void
uqTeuchosMatrixClass::cwSet(unsigned int rowId,unsigned int colId,const uqTeuchosMatrixClass& mat)
{
  UQ_FATAL_TEST_MACRO(rowId >= this->numRowsLocal(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::cwSet()",
                      "invalid rowId");

  UQ_FATAL_TEST_MACRO((rowId + mat.numRowsLocal()) > this->numRowsLocal(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::cwSet()",
                      "invalid vec.numRowsLocal()");

  UQ_FATAL_TEST_MACRO(colId >= this->numCols(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::cwSet()",
                      "invalid colId");

  UQ_FATAL_TEST_MACRO((colId + mat.numCols()) > this->numCols(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::cwSet()",
                      "invalid vec.numCols()");

  for (unsigned int i = 0; i < mat.numRowsLocal(); ++i) {
    for (unsigned int j = 0; j < mat.numCols(); ++j) {
      (*this)(rowId+i,colId+j) = mat(i,j);
    }
  }

  return;
}

// ---------------------------------------------------
//Kemelli tested 07/12/12
void
uqTeuchosMatrixClass::getColumn(unsigned int column_num, uqTeuchosVectorClass& column) const
{
  // Sanity checks
  UQ_FATAL_TEST_MACRO(column_num >= this->numCols(),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::getColumn",
                      "Specified row number not within range");

  UQ_FATAL_TEST_MACRO((column.sizeLocal() != this->numRowsLocal()),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::getColumn",
                      "column vector not same size as this matrix");

  // Temporary working pointer
  const double* temp_ptr ;
  
  // get the column_num- of matrix m_mat
  temp_ptr = m_mat[column_num];
  
  // Copy column from Teuchos matrix into our TeuchosVector object
  for (unsigned int i=0; i< column.sizeLocal();i++)
    column[i] = temp_ptr[i];
  
  return;
}

// ---------------------------------------------------
uqTeuchosVectorClass
uqTeuchosMatrixClass::getColumn(unsigned int column_num ) const
{
  // We rely on the void getColumn's to do sanity checks
  // So we don't do them here.
  uqTeuchosVectorClass column(m_env, m_map);
  this->getColumn( column_num, column );
  return column;
}


// ---------------------------------------------------
//Kemelli tested 07/12/12
void
uqTeuchosMatrixClass::setColumn(unsigned int column_num, const uqTeuchosVectorClass& column)
{
  this->resetLU();
  // Sanity checks
  UQ_FATAL_TEST_MACRO(column_num >= this->numCols(),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::setColumn",
                      "Specified column number not within range");

  UQ_FATAL_TEST_MACRO((column.sizeLocal() != this->numRowsLocal()),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::setColumn",
                      "column vector not same size as this matrix");
  
  for (unsigned int i =0; i < column.sizeLocal(); i++)
    m_mat(i,column_num) = column[i];
  
  return;
}

// ---------------------------------------------------
// tested 1/31/13
void
uqTeuchosMatrixClass::getRow(unsigned int row_num, uqTeuchosVectorClass& row) const
{
  // Sanity checks
  UQ_FATAL_TEST_MACRO(row_num >= this->numRowsLocal(),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::getRow",
                      "Specified row number not within range");

  UQ_FATAL_TEST_MACRO((row.sizeLocal() != this->numRowsLocal()),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::getRow",
                      "row vector not same size as this matrix");

  // Copy row from Teuchos matrix into our TeuchosVector object
  for (unsigned int i=0; i< row.sizeLocal();i++)
    row[i] = m_mat(row_num,i);
  
  return;
}

// ---------------------------------------------------
// tested 1/31/13
uqTeuchosVectorClass
uqTeuchosMatrixClass::getRow(unsigned int row_num ) const
{
  // We rely on the void getRow's to do sanity checks
  // So we don't do them here.
  uqTeuchosVectorClass row(m_env, m_map);
  this->getRow( row_num, row );
  return row;
}

// ---------------------------------------------------
// tested 1/31/13
void
uqTeuchosMatrixClass::setRow (const unsigned int row_num, const uqTeuchosVectorClass& row)
{
  // Sanity checks
  UQ_FATAL_TEST_MACRO(row_num >= this->numRowsLocal(),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::getRow",
                      "Specified row number not within range");

  UQ_FATAL_TEST_MACRO((row.sizeLocal() != this->numRowsLocal()),
                      env().fullRank(),
                      "uqTeuchosMatrixClass::getRow",
                      "row vector not same size as this matrix");

  // Copy our TeuchosVector object to our Teuchos Matrix
  for (unsigned int i=0; i< row.sizeLocal();i++)
     m_mat(row_num,i) = row[i] ;
  
  return;
}

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
void         
uqTeuchosMatrixClass::zeroLower(bool includeDiagonal)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((nRows != nCols),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::zeroLower()",
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

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
void 
uqTeuchosMatrixClass::zeroUpper(bool includeDiagonal)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();

  UQ_FATAL_TEST_MACRO((nRows != nCols),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::zeroUpper()",
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

// ---------------------------------------------------
void
uqTeuchosMatrixClass::filterSmallValues(double thresholdValue)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();
  
  for (unsigned int i = 0; i < nRows; ++i) {
    for (unsigned int j = 0; j < nCols; ++j) {
      double aux = (*this)(i,j);
      // If 'thresholdValue' is negative, no values will be filtered
      if ((aux < 0. ) && (-thresholdValue < aux)) {
        (*this)(i,j) = 0.;
      }      
      if ((aux > 0. ) && (thresholdValue > aux)) {
        (*this)(i,j) = 0.;
      }      
    }
  }
  return;
}

// ---------------------------------------------------
void
uqTeuchosMatrixClass::filterLargeValues(double thresholdValue)
{
  unsigned int nRows = this->numRowsLocal();
  unsigned int nCols = this->numCols();
  
  for (unsigned int i = 0; i < nRows; ++i) {
    for (unsigned int j = 0; j < nCols; ++j) {
      double aux = (*this)(i,j);
      // If 'thresholdValue' is negative, no values will be filtered
      if ( (aux < 0. ) && (-thresholdValue > aux)) 
          (*this)(i,j) = 0.;
      
      if ((aux > 0. ) && (thresholdValue < aux))
        (*this)(i,j) = 0.;
    }
  }
  return;
}


// ----------------------------------------------
// tested 1/31/13
void
uqTeuchosMatrixClass::fillWithTranspose(const uqTeuchosMatrixClass& mat)
{
  unsigned int nRows = mat.numRowsLocal();
  unsigned int nCols = mat.numCols();
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != nCols,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithTranspose()",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != nRows,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithTranspose()",
                      "inconsistent number of cols");

  for (unsigned int row = 0; row < nRows; ++row) {
    for (unsigned int col = 0; col < nCols; ++col) {
      (*this)(col,row) = mat(row,col);
    }
  }

  return;
}

// ---------------------------------------------------
// added 2/28/13
void
uqTeuchosMatrixClass::fillWithBlocksDiagonally(const std::vector<const uqTeuchosMatrixClass* >& matrices)
{
  unsigned int sumNumRowsLocals = 0;
  unsigned int sumNumCols       = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    sumNumRowsLocals += matrices[i]->numRowsLocal();
    sumNumCols       += matrices[i]->numCols();
  }
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != sumNumRowsLocals,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithBlocksDiagonally(const)",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != sumNumCols,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithBlocksDiagonally(const)",
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
// ---------------------------------------------------
// added 2/28/13
void
uqTeuchosMatrixClass::fillWithBlocksDiagonally(const std::vector<uqTeuchosMatrixClass* >& matrices)
{
  unsigned int sumNumRowsLocals = 0;
  unsigned int sumNumCols       = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    sumNumRowsLocals += matrices[i]->numRowsLocal();
    sumNumCols       += matrices[i]->numCols();
  }
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != sumNumRowsLocals,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithBlocksDiagonally()",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != sumNumCols,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithBlocksDiagonally()",
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
// ---------------------------------------------------
// added 2/28/13
void
uqTeuchosMatrixClass::fillWithBlocksHorizontally(const std::vector<const uqTeuchosMatrixClass* >& matrices)
{
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    UQ_FATAL_TEST_MACRO(this->numRowsLocal() != matrices[i]->numRowsLocal(),
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::fillWithBlocksHorizontally(const)",
                        "inconsistent local number of rows");
    sumNumCols += matrices[i]->numCols();
  }
  UQ_FATAL_TEST_MACRO(this->numCols() != sumNumCols,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithBlocksHorizontally(const)",
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

// ---------------------------------------------------
// added 2/28/13
void
uqTeuchosMatrixClass::fillWithBlocksHorizontally(const std::vector<uqTeuchosMatrixClass* >& matrices)
{
  unsigned int sumNumCols = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    UQ_FATAL_TEST_MACRO(this->numRowsLocal() != matrices[i]->numRowsLocal(),
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::fillWithBlocksHorizontally()",
                        "inconsistent local number of rows");
    sumNumCols += matrices[i]->numCols();
  }
  UQ_FATAL_TEST_MACRO(this->numCols() != sumNumCols,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithBlocksHorizontally()",
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

// ---------------------------------------------------
// added 2/28/13
void
uqTeuchosMatrixClass::fillWithBlocksVertically(const std::vector<const uqTeuchosMatrixClass* >& matrices)
{
  unsigned int sumNumRows = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    UQ_FATAL_TEST_MACRO(this->numCols() != matrices[i]->numCols(),
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::fillWithBlocksVertically(const)",
                        "inconsistent local number of cols");
    sumNumRows += matrices[i]->numRowsLocal();
  }
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != sumNumRows,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithBlocksVertically(const)",
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

// ---------------------------------------------------
// added 2/28/13
void
uqTeuchosMatrixClass::fillWithBlocksVertically(const std::vector<uqTeuchosMatrixClass* >& matrices)
{
  unsigned int sumNumRows = 0;
  for (unsigned int i = 0; i < matrices.size(); ++i) {
    UQ_FATAL_TEST_MACRO(this->numCols() != matrices[i]->numCols(),
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::fillWithBlocksVertically()",
                        "inconsistent local number of cols");
    sumNumRows += matrices[i]->numRowsLocal();
  }
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != sumNumRows,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillWithBlocksVertically()",
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

// ---------------------------------------------------
// added 2/28/13
void
uqTeuchosMatrixClass::fillWithTensorProduct(const uqTeuchosMatrixClass& mat1, const uqTeuchosMatrixClass& mat2)
{
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != (mat1.numRowsLocal() * mat2.numRowsLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillTensorProduct(mat and mat)",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != (mat1.numCols() * mat2.numCols()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillTensorProduct(mat and mat)",
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

// ---------------------------------------------------
// added 2/28/13
void
uqTeuchosMatrixClass::fillWithTensorProduct(const uqTeuchosMatrixClass& mat1, const uqTeuchosVectorClass& vec2)
{
  UQ_FATAL_TEST_MACRO(this->numRowsLocal() != (mat1.numRowsLocal() * vec2.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillTensorProduct(mat and vec)",
                      "inconsistent local number of rows");
  UQ_FATAL_TEST_MACRO(this->numCols() != (mat1.numCols() * 1),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::fillTensorProduct(mat and vec)",
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


// Miscellaneous methodos ----------------------------
// ---------------------------------------------------

void
uqTeuchosMatrixClass::mpiSum( const uqMpiCommClass& comm, uqTeuchosMatrixClass& M_global ) const
{
  // Sanity Checks
  UQ_FATAL_RC_MACRO(((this->numRowsLocal() != M_global.numRowsLocal()) ||
                     (this->numCols()      != M_global.numCols()     )),
		    env().fullRank(),
		    "uqTeuchosMatrixClass::mpiSum()",
		    "local and global matrices incompatible");

  /* TODO: Probably a better way to handle this unpacking/packing of data */
  int size = M_global.numRowsLocal()*M_global.numCols();
  std::vector<double> local( size, 0.0 );
  std::vector<double> global( size, 0.0 );

  int k;
  for( unsigned int i = 0; i < this->numRowsLocal(); i++ ) {
      for( unsigned int j = 0; j < this->numCols(); j++ ) {
	  k = i + j*M_global.numCols();
	  local[k] = (*this)(i,j);
	  }
  }

  comm.Allreduce((void*) &local[0], (void*) &global[0], size, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                 "uqTeuchosMatrixClass::mpiSum()",
                 "failed MPI.Allreduce()");

  for( unsigned int i = 0; i < this->numRowsLocal(); i++ ) {
      for( unsigned int j = 0; j < this->numCols(); j++ ) {
	  k = i + j*M_global.numCols();	  
	  M_global(i,j) = global[k];
	  }
  }

  return;
}

//--------------------------------------------------------
// tested 2/28/13
void
uqTeuchosMatrixClass::matlabLinearInterpExtrap(
  const uqTeuchosVectorClass& x1Vec,
  const uqTeuchosMatrixClass& y1Mat,
  const uqTeuchosVectorClass& x2Vec)
{
  UQ_FATAL_TEST_MACRO(x1Vec.sizeLocal() <= 1,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::matlabLinearInterpExtrap()",
                      "invalid 'x1' size");

  UQ_FATAL_TEST_MACRO(x1Vec.sizeLocal() != y1Mat.numRowsLocal(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::matlabLinearInterpExtrap()",
                      "invalid 'x1' and 'y1' sizes");

  UQ_FATAL_TEST_MACRO(x2Vec.sizeLocal() != this->numRowsLocal(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::matlabLinearInterpExtrap()",
                      "invalid 'x2' and 'this' sizes");

  UQ_FATAL_TEST_MACRO(y1Mat.numCols() != this->numCols(),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::matlabLinearInterpExtrap()",
                      "invalid 'y1' and 'this' sizes");

  uqTeuchosVectorClass y1Vec(x1Vec);
  uqTeuchosVectorClass y2Vec(x2Vec);
  for (unsigned int colId = 0; colId < y1Mat.numCols(); ++colId) {
    y1Mat.getColumn(colId,y1Vec);
    y2Vec.matlabLinearInterpExtrap(x1Vec,y1Vec,x2Vec);
    this->setColumn(colId,y2Vec);
  }

  return;
}

// I/O methodos --------------------------------------
// ---------------------------------------------------
// Kemelli 12/05/12 - tested
void
uqTeuchosMatrixClass::print(std::ostream& os) const
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
// ---------------------------------------------------

void
uqTeuchosMatrixClass::subReadContents(
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds)
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::subReadContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::subReadContents()",
                      "implemented just for sequential vectors for now");

  uqFilePtrSetStruct filePtrSet;
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

    // Read '=' sign
    *filePtrSet.ifsVar >> tmpString;

    UQ_FATAL_TEST_MACRO(tmpString != "=",
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::subReadContents()",
                        "string should be the '=' sign");

    // Read 'zeros(n_rows,n_cols)' string
    *filePtrSet.ifsVar >> tmpString;

    unsigned int posInTmpString = 6;

    // Isolate 'n_rows' in a string
    char nRowsString[tmpString.size()-posInTmpString+1];
    unsigned int posInRowsString = 0;
    do {
      UQ_FATAL_TEST_MACRO(posInTmpString >= tmpString.size(),
                          m_env.worldRank(),
                          "uqTeuchosMatrixClass::subReadContents()",
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
                          "uqTeuchosMatrixClass::subReadContents()",
                          "symbol ')' not found in first line of file");
      nColsString[posInColsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ')');
    nColsString[posInColsString] = '\0';

    // Convert 'n_rows' and 'n_cols' strings to numbers
    unsigned int numRowsInFile = (unsigned int) strtod(nRowsString,NULL);
    unsigned int numColsInFile = (unsigned int) strtod(nColsString,NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::subReadContents()"
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
                        "uqTeuchosMatrixClass::subReadContents()",
                        "size of vec in file is not big enough");

    // Check if [num of cols in file] == [num cols in current matrix]
    UQ_FATAL_TEST_MACRO(numColsInFile != nCols,
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::subReadContents()",
                        "number of parameters of vec in file is different than number of parameters in this vec object");

    // Code common to any core in a communicator
    unsigned int maxCharsPerLine = 64*nCols; // Up to about 60 characters to represent each parameter value

    unsigned int lineId = 0;
    while (lineId < idOfMyFirstLine) {
      filePtrSet.ifsVar->ignore(maxCharsPerLine,'\n');
      lineId++;
    };

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::subReadContents()"
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
                        "uqTeuchosMatrixClass::subReadContents()",
                        "in core 0, string should be the '=' sign");

    // Take into account the ' [' portion
    std::streampos tmpPos = filePtrSet.ifsVar->tellg();
    filePtrSet.ifsVar->seekg(tmpPos+(std::streampos)2);

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqTeuchosMatrixClass::subReadContents()"
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

// ---------------------------------------------------
void
uqTeuchosMatrixClass::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::subWriteContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::subWriteContents()",
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

//====================================================
// Private members
//====================================================

// ---------------------------------------------------
// Kemelli 12/05/12 - tested
void
uqTeuchosMatrixClass::copy(const uqTeuchosMatrixClass& src) //dummy
{
  this->resetLU();
  unsigned int i,j, nrows=src.numRowsLocal(), ncols=src.numCols();
  
  for(i=0; i< nrows ; i++)
    for (j = 0; j < ncols; j++)
      m_mat(i,j) = src(i,j);
    
  return;
}

// ---------------------------------------------------
void
uqTeuchosMatrixClass::resetLU()
{
  if (m_LU.numCols() >0 || m_LU.numRows() > 0) {
    m_LU.reshape(0,0); //Kemelli, 12/06/12, dummy    
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

  if (v_pivoting) {  //Kemelli added 12/09/12
    free(v_pivoting);
    v_pivoting = NULL;
    
  }
  m_signum = 0;
  m_isSingular = false;

  return;
}

// ---------------------------------------------------
// multiply this matrix by vector x and store in vector y
// checked 12/10/12
void
uqTeuchosMatrixClass::multiply(const uqTeuchosVectorClass& x, uqTeuchosVectorClass& y) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != x.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::multiply(), vector return void",
                      "matrix and x have incompatible sizes");

  UQ_FATAL_TEST_MACRO((this->numRowsLocal() != y.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosMatrixClass::multiply(), vector return void",
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

// ---------------------------------------------------
// Implemented(finally) and checked 1/10/13
int
uqTeuchosMatrixClass::internalSvd() const
{
  if (m_svdColMap == NULL) {
    int nRows = (int) this->numRowsLocal();
    int nCols = (int) this->numCols();
    UQ_FATAL_TEST_MACRO(nRows < nCols,
                        m_env.worldRank(),
                        "uqTeuchosMatrixClass::internalSvd()",
                        "LAPACK/Teuchos only supports cases where nRows >= nCols");

    m_svdColMap = new uqMapClass(this->numCols(),0,this->map().Comm()); // see 'uqVectorSpaceClass<.,.>::newMap()'
  //in src/basic/src/uqTeuchosVectorSpace.C //old comment already existent in uqGslMatrixClass
    m_svdUmat   = new uqTeuchosMatrixClass(*this); // Yes, 'this'
    m_svdSvec   = new uqTeuchosVectorClass(m_env,*m_svdColMap);
    m_svdVmat   = new uqTeuchosMatrixClass(*m_svdSvec);
    m_svdVTmat  = new uqTeuchosMatrixClass(*m_svdSvec);
    
    int minRowsCols, maxRowsCols;
    
    if (nRows>=nCols) { minRowsCols = nCols; maxRowsCols = nRows; } else { minRowsCols = nRows; maxRowsCols = nCols; }

    char jobu, jobvt;
    int  lwork, info;
    double  *work, *rwork;
  
    jobu = 'S';
    jobvt = 'S';
  
    lwork = 15*maxRowsCols; // Set up the work array, larger than needed.
    work = new double[lwork];

    int aux1= 5*minRowsCols+7, aux2= 2*maxRowsCols+2*minRowsCols+1;
    int aux_dim;
  
    if (aux1>=aux2) { aux_dim = minRowsCols*aux1; } else {aux_dim = minRowsCols*aux2; }
  
    rwork = new double[aux_dim];
  
    Teuchos::LAPACK<int, double> lapack;
    
    lapack.GESVD(jobu,jobvt,m_mat.numRows(),m_mat.numCols(),m_mat.values(),m_mat.stride(),
	        m_svdSvec->values(),m_svdUmat->values(),m_svdUmat->stride(),m_svdVTmat->values(),
		m_svdVTmat->stride(),&work[0],lwork,&rwork[0],&info);
  }
  return 0;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++
// Operators outside class definition
//++++++++++++++++++++++++++++++++++++++++++++++++++++

uqTeuchosMatrixClass operator*(double a, const uqTeuchosMatrixClass& mat)
{
  uqTeuchosMatrixClass answer(mat);
  answer *= a;
  return answer;
}

// ---------------------------------------------------
uqTeuchosVectorClass operator*(const uqTeuchosMatrixClass& mat, const uqTeuchosVectorClass& vec)
{
  return mat.multiply(vec);
}

// ---------------------------------------------------
uqTeuchosMatrixClass operator*(const uqTeuchosMatrixClass& m1, const uqTeuchosMatrixClass& m2)
{
  unsigned int m1Rows = m1.numRowsLocal();
  unsigned int m1Cols = m1.numCols();
  unsigned int m2Rows = m2.numRowsLocal();
  unsigned int m2Cols = m2.numCols();

  UQ_FATAL_TEST_MACRO((m1Cols != m2Rows),
                      m1.env().worldRank(),
                      "uqTeuchosMatrixClass operator*(matrix,matrix)",
                      "different sizes m1Cols and m2Rows");

  uqTeuchosMatrixClass mat(m1.env(),m1.map(),m2Cols);
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

// ---------------------------------------------------
uqTeuchosMatrixClass operator+(const uqTeuchosMatrixClass& m1, const uqTeuchosMatrixClass& m2)
{
  uqTeuchosMatrixClass answer(m1);
  answer += m2;
  return answer;
}

// ---------------------------------------------------
uqTeuchosMatrixClass operator-(const uqTeuchosMatrixClass& m1, const uqTeuchosMatrixClass& m2)
{
  uqTeuchosMatrixClass answer(m1);
  answer -= m2;
  return answer;
}

// ---------------------------------------------------
uqTeuchosMatrixClass matrixProduct(const uqTeuchosVectorClass& v1, const uqTeuchosVectorClass& v2)
{
  unsigned int nRows = v1.sizeLocal();
  unsigned int nCols = v2.sizeLocal();
  uqTeuchosMatrixClass answer(v1.env(),v1.map(),nCols);

  for (unsigned int i = 0; i < nRows; ++i) {
    double value1 = v1[i];
    for (unsigned int j = 0; j < nCols; ++j) {
      answer(i,j) = value1*v2[j];
    }
  }

  return answer;
}

// ---------------------------------------------------
uqTeuchosMatrixClass leftDiagScaling(const uqTeuchosVectorClass& vec, const uqTeuchosMatrixClass& mat)
{
  unsigned int vSize = vec.sizeLocal();
  unsigned int mRows = mat.numRowsLocal();
  unsigned int mCols = mat.numCols();

  UQ_FATAL_TEST_MACRO((vSize != mRows),
                      mat.env().worldRank(),
                      "uqTeuchosMatrixClass leftDiagScaling(vector,matrix)",
                      "size of vector is different from the number of rows in matrix");

  UQ_FATAL_TEST_MACRO((mCols != mRows),
                      mat.env().worldRank(),
                      "uqTeuchosMatrixClass leftDiagScaling(vector,matrix)",
                      "routine currently works for square matrices only");

  uqTeuchosMatrixClass answer(mat);
  for (unsigned int i = 0; i < mRows; ++i) {
    double vecValue = vec[i];
    for (unsigned int j = 0; j < mCols; ++j) {
      answer(i,j) *= vecValue;
    }
  }

  return answer;
}

// ---------------------------------------------------
uqTeuchosMatrixClass rightDiagScaling(const uqTeuchosMatrixClass& mat, const uqTeuchosVectorClass& vec)
{
  unsigned int vSize = vec.sizeLocal();
  unsigned int mRows = mat.numRowsLocal();
  unsigned int mCols = mat.numCols();

  UQ_FATAL_TEST_MACRO((vSize != mCols),
                      mat.env().worldRank(),
                      "uqTeuchosMatrixClass rightDiagScaling(matrix,vector)",
                      "size of vector is different from the number of cols in matrix");

  UQ_FATAL_TEST_MACRO((mCols != mRows),
                      mat.env().worldRank(),
                      "uqTeuchosMatrixClass rightDiagScaling(matrix,vector)",
                      "routine currently works for square matrices only");

  uqTeuchosMatrixClass answer(mat);
  for (unsigned int j = 0; j < mCols; ++j) {
    double vecValue = vec[j];
    for (unsigned int i = 0; i < mRows; ++i) {
      answer(i,j) *= vecValue;
    }
  }

  return answer;
}

// ---------------------------------------------------
std::ostream&
operator<<(std::ostream& os, const uqTeuchosMatrixClass& obj)
{
  obj.print(os);

  return os;
}

}  // End namespace QUESO

#endif // ifdef QUESO_HAS_TRILINOS
