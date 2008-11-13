/* libs/basic/src/uqTrilinosVector.C
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

#include <uqTrilinosVector.h>
#include <Epetra_MpiComm.h>
#include <uqDefines.h>

uqTrilinosVectorClass::uqTrilinosVectorClass()
  :
  uqVectorClass(),
  m_map        (*(new Epetra_Map( 1,0,*(new Epetra_MpiComm(MPI_COMM_WORLD)) ) ))
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqTrilinosVectorClass::constructor(), default",
                      "should not be used by user");
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map)
  :
  uqVectorClass(env),
  m_map        (map),
  //m_vec      (new Epetra_Vector(map,true))
  m_vec        (new Epetra_SerialDenseMatrix(map.NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.rank(),
                      "uqTrilinosVectorClass::constructor()",
                      "null vector generated");
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map, double d1, double d2, unsigned int size)
  :
  uqVectorClass(env),
  m_map        (map),
  //m_vec      (new Epetra_Vector(map,true))
  m_vec        (new Epetra_SerialDenseMatrix(map.NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.rank(),
                      "uqTrilinosVectorClass::constructor(), linspace",
                      "null vector generated");

  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::constructor(), linspace",
                    "failed");
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqTrilinosVectorClass& v, double d1, double d2, unsigned int size)
  :
  uqVectorClass(v.env()),
  m_map        (v.map()),
  //m_vec      (new Epetra_Vector(v.map(),true))
  m_vec        (new Epetra_SerialDenseMatrix(v.map().NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.rank(),
                      "uqTrilinosVectorClass::constructor(), linspace",
                      "null vector generated");

  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::constructor(), linspace",
                    "failed");
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqTrilinosVectorClass& v)
  :
  uqVectorClass(v.env()),
  m_map        (v.map()),
  //m_vec      (new Epetra_Vector(v.map(),true))
  m_vec        (new Epetra_SerialDenseMatrix(v.map().NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.rank(),
                      "uqTrilinosVectorClass::constructor(), copy",
                      "null vector generated");
  this->uqVectorClass::copy(v);
  this->copy(v);
}

uqTrilinosVectorClass::~uqTrilinosVectorClass()
{
  if (m_vec) delete m_vec;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator=(const uqTrilinosVectorClass& obj)
{
  this->copy(obj);
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator*=(double a)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::operator*=()",
                    "failed");
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator/=(double a)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::operator/=()",
                    "failed");
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator*=(const uqTrilinosVectorClass& rhs)
{
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator/=(const uqTrilinosVectorClass& rhs)
{
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator+=(const uqTrilinosVectorClass& rhs)
{
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator-=(const uqTrilinosVectorClass& rhs)
{
  return *this;
}

double&
uqTrilinosVectorClass::operator[](unsigned int i)
{
  return (*m_vec)(i,0);
}

const double&
uqTrilinosVectorClass::operator[](unsigned int i) const
{
  return (*m_vec)(i,0);
}

void
uqTrilinosVectorClass::copy(const uqTrilinosVectorClass& src)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::copy()",
                    "failed");

  return;
}

unsigned int
uqTrilinosVectorClass::size() const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::size()",
                    "failed");

  return 0;
}

double
uqTrilinosVectorClass::norm2Sq() const
{
  return 0.;
}

double
uqTrilinosVectorClass::norm2() const
{
  return 0.;
}

double
uqTrilinosVectorClass::sumOfComponents() const
{
  double result = 0.;

  return result;
}

void
uqTrilinosVectorClass::cwSet(double value)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::set()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::cwSetGaussian(gsl_rng* rng, double mean, double stdDev)
{
  return;
}

void
uqTrilinosVectorClass::cwInvert()
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::invert()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::sort()
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::sort()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::print(std::ostream& os) const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.rank(),
                    "uqTrilinosVectorClass::print()",
                    "failed");
  return;
}

//Epetra_Vector*
Epetra_SerialDenseMatrix*
uqTrilinosVectorClass::data() const
{
  return m_vec;
}

bool
uqTrilinosVectorClass::atLeastOneComponentSmallerThan(const uqTrilinosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->size() != rhs.size()),
                      m_env.rank(),
                      "uqTrilinosVectorClass::atLeastOneComponentSmallerThan()",
                      "vectors have different sizes");

  bool result = false;

  return result;
}

bool
uqTrilinosVectorClass::atLeastOneComponentBiggerThan (const uqTrilinosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->size() != rhs.size()),
                      m_env.rank(),
                      "uqTrilinosVectorClass::atLeastOneComponentBiggerThan()",
                      "vectors have different sizes");

  bool result = false;

  return result;
}

std::ostream&
operator<<(std::ostream& os, const uqTrilinosVectorClass& obj)
{
  obj.print(os);

  return os;
}

uqTrilinosVectorClass operator/(const double a, const uqTrilinosVectorClass& x)
{
  uqTrilinosVectorClass answer(x);
  answer.cwInvert();
  answer *= a;

  return answer;
}

uqTrilinosVectorClass operator/(const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y)
{
  uqTrilinosVectorClass answer(x);
  answer /= y;

  return answer;
}

uqTrilinosVectorClass operator*(const double a, const uqTrilinosVectorClass& x)
{
  uqTrilinosVectorClass answer(x);
  answer *= a;

  return answer;
}

uqTrilinosVectorClass operator*(const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y)
{
  uqTrilinosVectorClass answer(x);
  answer *= y;

  return answer;
}

double scalarProduct(const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y)
{
  unsigned int size1 = x.size();
  unsigned int size2 = y.size();
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      x.env().rank(),
                      "scalarProduct()",
                      "different sizes of x and y");

  double result = 0.;
  //for (unsigned int i = 0; i < size1; ++i) {
  //  result += x[i]*y[i];
  //}

  return result;
}

uqTrilinosVectorClass operator+(const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y)
{
  uqTrilinosVectorClass answer(x);
  answer += y;

  return answer;
}

uqTrilinosVectorClass operator-(const uqTrilinosVectorClass& x, const uqTrilinosVectorClass& y)
{
  uqTrilinosVectorClass answer(x);
  answer -= y;

  return answer;
}

const Epetra_Map&
uqTrilinosVectorClass::map() const
{
  return m_map;
  //return (const Epetra_Map&) (m_vec->Map());
}
#if 0
int
uqTrilinosVectorClass::rank() const
{
  return m_vec->Map().Comm().MyPID();
}
#endif
