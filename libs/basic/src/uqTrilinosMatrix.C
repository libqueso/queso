/* libs/basic/src/uqTrilinosMatrix.C
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

#include <uqTrilinosMatrix.h>
#include <uqTrilinosVector.h>
#include <Epetra_MpiComm.h>
#include <uqDefines.h>

uqTrilinosMatrixClass::uqTrilinosMatrixClass()
  :
  uqMatrixClass(),
  m_map        (*(new Epetra_Map( 1,0,*(new Epetra_MpiComm(MPI_COMM_WORLD)) ) ))
{
}
 
uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqBaseEnvironmentClass& env,
  const Epetra_Map&         map,
  unsigned int              numCols)
  :
  uqMatrixClass(env),
  m_map        (map),
  //m_mat(new Epetra_CrsMatrix(Copy,map,numCols))
  m_mat(new Epetra_SerialDenseMatrix(map.NumGlobalElements(),numCols))
{
  UQ_FATAL_TEST_MACRO((m_mat == NULL),
                      m_env.rank(),
                      "uqTrilinosMatrixClass::constructor()",
                      "null matrix generated");
}
 
uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqBaseEnvironmentClass& env,
  const Epetra_Map&         map,
  double                    diagValue)
  :
  uqMatrixClass(env),
  m_map        (map)
{
}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(
  const uqTrilinosVectorClass& v,
  double                       diagValue)
  :
  uqMatrixClass(v.env()),
  m_map        (v.map())
{
}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(const uqTrilinosVectorClass& v)
  :
  uqMatrixClass(v.env()),
  m_map        (v.map())
{
}

uqTrilinosMatrixClass::uqTrilinosMatrixClass(const uqTrilinosMatrixClass& B)
  :
  uqMatrixClass(B.env()),
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
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator/=(double a)
{
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator+=(const uqTrilinosMatrixClass& rhs)
{
  return *this;
}

uqTrilinosMatrixClass&
uqTrilinosMatrixClass::operator-=(const uqTrilinosMatrixClass& rhs)
{
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
  return;
}

unsigned int
uqTrilinosMatrixClass::numRows() const
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
  return;
}

void
uqTrilinosMatrixClass::zeroUpper(bool includeDiagonal)
{
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
  UQ_FATAL_TEST_MACRO((this->numCols() != x.size()),
                      m_env.rank(),
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
  UQ_FATAL_TEST_MACRO((this->numCols() != b.size()),
                      m_env.rank(),
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
  UQ_FATAL_TEST_MACRO((this->numCols() != b.size()),
                      m_env.rank(),
                      "uqTrilinosMatrixClass::multiply(), return void",
                      "matrix and rhs have incompatible sizes");

  UQ_FATAL_TEST_MACRO((x.size() != b.size()),
                      m_env.rank(),
                      "uqTrilinosMatrixClass::multiply(), return void",
                      "solution and rhs have incompatible sizes");

  return;
}

void
uqTrilinosMatrixClass::print(std::ostream& os) const
{
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
  unsigned int m1Rows = m1.numRows();
  unsigned int m1Cols = m1.numCols();
  unsigned int m2Rows = m2.numRows();
  unsigned int m2Cols = m2.numCols();

  UQ_FATAL_TEST_MACRO((m1Cols != m2Rows),
                      m1.env().rank(),
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
  unsigned int numRows = v1.size();
  unsigned int numCols = v2.size();
  uqTrilinosMatrixClass answer(v1.env(),v1.map(),numRows,numCols);

  //for (unsigned int i = 0; i < numRows; ++i) {
  //  double value1 = v1[i];
  //  for (unsigned int j = 0; j < numCols; ++j) {
  //    answer(i,j) = value1*v2[j];
  //  }
  //}

  return answer;
}

uqTrilinosMatrixClass diagScaling(const uqTrilinosVectorClass& vec, const uqTrilinosMatrixClass& mat)
{
  unsigned int vSize = vec.size();
  unsigned int mRows = mat.numRows();
  unsigned int mCols = mat.numCols();

  UQ_FATAL_TEST_MACRO((vSize != mRows),
                      mat.env().rank(),
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
uqTrilinosMatrixClass::rank() const
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
