//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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

#include <uqTrilinosMatrix.h>
#include <Epetra_MpiComm.h>
#include <uqDefines.h>

uqTrilinosMatrixClass::uqTrilinosMatrixClass()
  :
  uqMatrixClass(),
  m_map        (*(new Epetra_Map( 1,0,*(new Epetra_MpiComm(MPI_COMM_WORLD)) ) ))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqTrilinosMatrixClass::constructor(), default",
                      "should not be used by user");
}
 
uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqBaseEnvironmentClass& env,
  const Epetra_Map&         map,
  unsigned int              numCols)
  :
  uqMatrixClass(env, map),
  m_map        (map),
  //m_mat(new Epetra_CrsMatrix(Copy,map,numCols))
  m_mat(new Epetra_SerialDenseMatrix(map.NumGlobalElements(),numCols))
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.fullRank(),
                      "uqTrilinosMatrixClass::constructor()",
                      "null matrix generated");
}
 
uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqBaseEnvironmentClass& env,
  const Epetra_Map&         map,
  unsigned int              numRowsLocal,
  unsigned int              numCols)
  :
  uqMatrixClass(env, map),
  m_map        (map),
  m_mat(new Epetra_SerialDenseMatrix(numRowsLocal,numCols))
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.fullRank(),
                      "uqTrilinosMatrixClass::constructor()",
                      "null matrix generated");
}
 
uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqBaseEnvironmentClass& env,
  const Epetra_Map&         map,
  double                    diagValue)
  :
  uqMatrixClass(env, map),
  m_map        (map)
{
  double x = diagValue; x += 1.; // just to avoid icpc warnings
}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqTrilinosVectorClass& v,
  double                       diagValue)
  :
  uqMatrixClass(v.env(), v.map()),
  m_map        (v.map())
{
  double x = diagValue; x += 1.; // just to avoid icpc warnings
}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(const uqTrilinosVectorClass& v)
  :
  uqMatrixClass(v.env(), v.map()),
  m_map        (v.map())
{
}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(const uqTrilinosMatrixClass& B)
  :
  uqMatrixClass(B.env(), B.map()),
  m_map        (B.map())
{
}

uqTrilinosMatrixClass::~uqTrilinosMatrixClass()
{
  if (m_mat) delete m_mat;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator=(const uqTrilinosMatrixClass& obj)
{
  this->copy(obj);
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator*=(double a)
{
  double x = a; x += 1.; // just to avoid icpc warnings
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator/=(double a)
{
  double x = a; x += 1.; // just to avoid icpc warnings
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator+=(const uqTrilinosMatrixClass& rhs)
{
  double tmpA = rhs(0,0); tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator-=(const uqTrilinosMatrixClass& rhs)
{
  double tmpA = rhs(0,0); tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

double&
uqTrilinosMatrixClass::operator()(unsigned int i, unsigned int j)
{
  return (*m_mat)(i,j);
}

const double&
uqTrilinosMatrixClass::operator()(unsigned int i, unsigned int j) const
{
  return (*m_mat)(i,j);
}

void
uqTrilinosMatrixClass::copy(const uqTrilinosMatrixClass& src)
{
  double tmpA = src(0,0); tmpA += 1.; // Just to avoid icpc warnings
  return;
}

unsigned int
uqTrilinosMatrixClass::numRowsLocal() const
{
  return 0;
}

unsigned int
uqTrilinosMatrixClass::numRowsGlobal() const
{
  return 0;
}

unsigned int
uqTrilinosMatrixClass::numCols() const
{
  return 0;
}

int
uqTrilinosMatrixClass::chol()
{
  return 0;
}

void
uqTrilinosMatrixClass::zeroLower(bool includeDiagonal)
{
  double tmpA = (double) includeDiagonal; tmpA += 1.; // Just to avoid icpc warnings
  return;
}

void
uqTrilinosMatrixClass::zeroUpper(bool includeDiagonal)
{
  double tmpA = (double) includeDiagonal; tmpA += 1.; // Just to avoid icpc warnings
  return;
}

uqTrilinosMatrixClass
uqTrilinosMatrixClass::transpose() const
{
  return *this;
}

uqTrilinosVectorClass
uqTrilinosMatrixClass::multiply(
  const uqTrilinosVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != x.sizeLocal()),
                      m_env.fullRank(),
                      "uqTrilinosMatrixClass::multiply(), vector",
                      "matrix and x have incompatible sizes");

  uqTrilinosVectorClass y(this->env(),this->map());
  m_mat->Multiply( false,*(x.data()),*(y.data()) );

  return y;
}

uqTrilinosVectorClass
uqTrilinosMatrixClass::invertMultiply(
  const uqTrilinosVectorClass& b) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.fullRank(),
                      "uqTrilinosMatrixClass::invertMultiply(), return vector",
                      "matrix and rhs have incompatible sizes");

  uqTrilinosVectorClass x(m_env,this->map());
  this->invertMultiply(b,x);

  return x;
}

void
uqTrilinosMatrixClass::invertMultiply(
  const uqTrilinosVectorClass& b,
        uqTrilinosVectorClass& x) const
{
  UQ_FATAL_TEST_MACRO((this->numCols() != b.sizeLocal()),
                      m_env.fullRank(),
                      "uqTrilinosMatrixClass::multiply(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.sizeLocal() != b.sizeLocal()),
                      m_env.fullRank(),
                      "uqTrilinosMatrixClass::multiply(), return void",
                      "solution and rhs have incompatible sizes");

  return;
}

void
uqTrilinosMatrixClass::print(std::ostream& os) const
{
  os.flush(); // just to avoid icpc warnings
  return;
}

std::ostream&
operator<<(std::ostream& os, const uqTrilinosMatrixClass& obj)
{
  obj.print(os);

  return os;
}

uqTrilinosMatrixClass operator*(double a, const uqTrilinosMatrixClass& mat)
{
  uqTrilinosMatrixClass answer(mat);
  answer *= a;
  return answer;
}

uqTrilinosVectorClass operator*(const uqTrilinosMatrixClass& mat, const uqTrilinosVectorClass& vec)
{
  return mat.multiply(vec);
}

uqTrilinosMatrixClass operator*(const uqTrilinosMatrixClass& m1, const uqTrilinosMatrixClass& m2)
{
  unsigned int m1Rows = m1.numRowsLocal();
  unsigned int m1Cols = m1.numCols();
  unsigned int m2Rows = m2.numRowsLocal();
  unsigned int m2Cols = m2.numCols();

  UQ_FATAL_TEST_MACRO((m1Cols != m2Rows),
                      m1.env().fullRank(),
                      "uqTrilinosMatrixClass operator*(matrix,matrix)",
                      "different sizes m1Cols and m2Rows");

  uqTrilinosMatrixClass mat(m1.env(),m1.map(),m1Rows,m2Cols);

  //unsigned int commonSize = m1Cols;
  //for (unsigned int row1 = 0; row1 < m1Rows; ++row1) {
  //  for (unsigned int col2 = 0; col2 < m2Cols; ++col2) {
  //    double result = 0.;
  //    for (unsigned int k = 0; k < commonSize; ++k) {
  //      result += m1(row1,k)*m2(k,col2);
  //    }
  //    mat(row1,col2) = result;
  //  }
  //}

  return mat;
}

uqTrilinosMatrixClass operator+(const uqTrilinosMatrixClass& m1, const uqTrilinosMatrixClass& m2)
{
  uqTrilinosMatrixClass answer(m1);
  answer += m2;
  return answer;
}

uqTrilinosMatrixClass matrixProduct(const uqTrilinosVectorClass& v1, const uqTrilinosVectorClass& v2)
{
  unsigned int numRowsLocal = v1.sizeLocal();
  unsigned int numCols = v2.sizeLocal();
  uqTrilinosMatrixClass answer(v1.env(),v1.map(),numRowsLocal,numCols);

  //for (unsigned int i = 0; i < numRowsLocal; ++i) {
  //  double value1 = v1[i];
  //  for (unsigned int j = 0; j < numCols; ++j) {
  //    answer(i,j) = value1*v2[j];
  //  }
  //}

  return answer;
}

uqTrilinosMatrixClass diagScaling(const uqTrilinosVectorClass& vec, const uqTrilinosMatrixClass& mat)
{
  unsigned int vSize = vec.sizeLocal();
  unsigned int mRows = mat.numRowsLocal();
  //unsigned int mCols = mat.numCols();

  UQ_FATAL_TEST_MACRO((vSize != mRows),
                      mat.env().fullRank(),
                      "uqTrilinosMatrixClass diagScaling(vector,matrix)",
                      "size of vector is different from the number of rows in matrix");

  uqTrilinosMatrixClass answer(mat);
  //for (unsigned int i = 0; i < mRows; ++i) {
  //  double vecValue = vec[i];
  //  for (unsigned int j = 0; j < mCols; ++j) {
  //    answer(i,j) *= vecValue;
  //  }
  //}

  return answer;
}

#if 0
int
uqTrilinosMatrixClass::fullRank() const
{
  return this->map().Comm().MyPID();
}

double
uqTrilinosMatrixClass::get(unsigned int i, unsigned int j) const
{
  return 0.;
}

void
uqTrilinosMatrixClass::set(unsigned int i, unsigned int j, double value)
{
  return;
}

void
uqTrilinosMatrixClass::scale(double alpha)
{
  return;
}

void
uqTrilinosMatrixClass::add(const uqTrilinosMatrixClass* B)
{
  return;
}

void
uqTrilinosMatrixClass::sub(const uqTrilinosMatrixClass* B)
{
  return;
}
#endif
