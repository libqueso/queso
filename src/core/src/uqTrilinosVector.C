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

#include <uqTrilinosVector.h>
#ifdef QUESO_HAS_TRILINOS

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>

uqTrilinosVectorClass::uqTrilinosVectorClass()
  :
  uqVectorClass(),
  m_map        (*(new uqMapClass( 1,0,*(new uqMpiCommClass(*(new uqEmptyEnvironmentClass()),uqRawValue_MPI_COMM_SELF)) ) )) // avoid using MPI_COMM_WORLD
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqTrilinosVectorClass::constructor(), default",
                      "should not be used by user");
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map)
  :
  uqVectorClass(env, map),
  m_map        (map),
  //m_vec      (new Epetra_Vector(map,true))
  m_vec        (new Epetra_SerialDenseMatrix(map.NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqTrilinosVectorClass::constructor()",
                      "null vector generated");
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double value)
  :
  uqVectorClass(env,map),
  m_map        (map),
  //m_vec      (new Epetra_Vector(map,true))
  m_vec        (new Epetra_SerialDenseMatrix(map.NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqTrilinosVectorClass::constructor()",
                      "null vector generated");
  this->cwSet(value);
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double d1, double d2, unsigned int size)
  :
  uqVectorClass(env, map),
  m_map        (map),
  //m_vec      (new Epetra_Vector(map,true))
  m_vec        (new Epetra_SerialDenseMatrix(map.NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqTrilinosVectorClass::constructor(), linspace",
                      "null vector generated");

  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::constructor(), linspace",
                    "failed");

  double x = d1+d2;   // just to avoid icpc warnings
  x += (double) size; // just to avoid icpc warnings
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqTrilinosVectorClass& v, double d1, double d2, unsigned int size)
  :
  uqVectorClass(v.env(), v.map()),
  m_map        (v.map()),
  //m_vec      (new Epetra_Vector(v.map(),true))
  m_vec        (new Epetra_SerialDenseMatrix(v.map().NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqTrilinosVectorClass::constructor(), linspace",
                      "null vector generated");

  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::constructor(), linspace",
                    "failed");

  double x = d1+d2;   // just to avoid icpc warnings
  x += (double) size; // just to avoid icpc warnings
}

uqTrilinosVectorClass::uqTrilinosVectorClass(const uqTrilinosVectorClass& v)
  :
  uqVectorClass(v.env(), v.map()),
  m_map        (v.map()),
  //m_vec      (new Epetra_Vector(v.map(),true))
  m_vec        (new Epetra_SerialDenseMatrix(v.map().NumGlobalElements(),1))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
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
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::operator*=()",
                    "failed");
  double tmpA = a; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator/=(double a)
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::operator/=()",
                    "failed");
  double tmpA = a; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator*=(const uqTrilinosVectorClass& rhs)
{
  double tmpA = rhs[0]; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator/=(const uqTrilinosVectorClass& rhs)
{
  double tmpA = rhs[0]; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator+=(const uqTrilinosVectorClass& rhs)
{
  double tmpA = rhs[0]; tmpA += 1.; // Just to avoid icpc warnings
  return *this;
}

uqTrilinosVectorClass&
uqTrilinosVectorClass::operator-=(const uqTrilinosVectorClass& rhs)
{
  double tmpA = rhs[0]; tmpA += 1.; // Just to avoid icpc warnings
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
                    src.env().fullRank(),
                    "uqTrilinosVectorClass::copy()",
                    "failed");

  return;
}

unsigned int
uqTrilinosVectorClass::sizeLocal() const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::sizeLocal()",
                    "failed");

  return 0;
}

unsigned int
uqTrilinosVectorClass::sizeGlobal() const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::sizeGlobal()",
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
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::set()",
                    "failed");
  double tmpA = value; tmpA += 1.; // Just to avoid icpc warnings
  return;
}

void
uqTrilinosVectorClass::cwSetGaussian(const gsl_rng* rng, double mean, double stdDev)
{
  double tmpA = mean + stdDev;          // Just to avoid icpc warnings
  tmpA += gsl_ran_gaussian(rng,stdDev); // Just to avoid icpc warnings
  return;
}

void
uqTrilinosVectorClass::cwInvert()
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::invert()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::sort()
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::sort()",
                    "failed");
  return;
}

void
uqTrilinosVectorClass::print(std::ostream& os) const
{
  UQ_FATAL_RC_MACRO(UQ_INCOMPLETE_IMPLEMENTATION_RC,
                    m_env.fullRank(),
                    "uqTrilinosVectorClass::print()",
                    "failed");
  os.flush(); // just to avoid icpc warnings
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
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.fullRank(),
                      "uqTrilinosVectorClass::atLeastOneComponentSmallerThan()",
                      "vectors have different sizes");

  bool result = false;

  return result;
}

bool
uqTrilinosVectorClass::atLeastOneComponentBiggerThan (const uqTrilinosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.fullRank(),
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

uqTrilinosVectorClass operator/(double a, const uqTrilinosVectorClass& x)
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

uqTrilinosVectorClass operator*(double a, const uqTrilinosVectorClass& x)
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
  unsigned int size1 = x.sizeLocal();
  unsigned int size2 = y.sizeLocal();
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      x.env().fullRank(),
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

const uqMapClass&
uqTrilinosVectorClass::map() const
{
  return m_map;
  //return (const uqMapClass&) (m_vec->Map());
}
#if 0
int
uqTrilinosVectorClass::fullRank() const
{
  return m_vec->Map().Comm().MyPID();
}
#endif

#endif // #ifdef QUESO_HAS_TRILINOS
