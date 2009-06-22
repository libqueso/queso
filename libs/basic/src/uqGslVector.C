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

#include <uqGslVector.h>
#include <uqDefines.h>
#include <gsl/gsl_sort_vector.h>
#include <math.h>

uqGslVectorClass::uqGslVectorClass()
  :
  uqVectorClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.fullRank(),
                      "uqGslVectorClass::constructor(), default",
                      "should not be used by user");
}

uqGslVectorClass::uqGslVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map)
  :
  uqVectorClass(env,map),
  m_vec        (gsl_vector_calloc(map.NumGlobalElements()))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqGslVectorClass::constructor()",
                      "null vector generated");
}

uqGslVectorClass::uqGslVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map, double value)
  :
  uqVectorClass(env,map),
  m_vec(gsl_vector_calloc(map.NumGlobalElements()))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqGslVectorClass::constructor()",
                      "null vector generated");
  this->cwSet(value);
}

uqGslVectorClass::uqGslVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const Epetra_Map& map)
  :
  uqVectorClass(env,map),
  m_vec(gsl_vector_calloc(map.NumGlobalElements()))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqGslVectorClass::constructor(), linspace",
                      "null vector generated");

  for (unsigned int i = 0; i < m_vec->size; ++i) {
    double alpha = (double) i / ((double) m_vec->size - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }
}

uqGslVectorClass::uqGslVectorClass(const uqGslVectorClass& v, double d1, double d2)
  :
  uqVectorClass(v.env(),v.map()),
  m_vec(gsl_vector_calloc(v.sizeLocal()))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqGslVectorClass::constructor(), linspace",
                      "null vector generated");

  for (unsigned int i = 0; i < m_vec->size; ++i) {
    double alpha = (double) i / ((double) m_vec->size - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }
}

uqGslVectorClass::uqGslVectorClass(const uqGslVectorClass& v)
  :
  uqVectorClass(v.env(),v.map()),
  m_vec(gsl_vector_calloc(v.sizeLocal()))
{
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.fullRank(),
                      "uqGslVectorClass::constructor(), copy",
                      "null vector generated");
  this->uqVectorClass::copy(v);
  this->copy(v);
}

uqGslVectorClass::~uqGslVectorClass()
{
  if (m_vec) gsl_vector_free(m_vec);
}

unsigned int
uqGslVectorClass::numberOfProcessorsRequiredForStorage() const
{
  return 1;
}

uqGslVectorClass&
uqGslVectorClass::operator=(const uqGslVectorClass& rhs)
{
  this->copy(rhs);
  return *this;
}

uqGslVectorClass&
uqGslVectorClass::operator*=(double a)
{
  int iRC;
  iRC = gsl_vector_scale(m_vec,a);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslVectorClass::operator*=()",
                    "failed");
  return *this;
}

uqGslVectorClass&
uqGslVectorClass::operator/=(double a)
{
  *this *= (1./a);

  return *this;
}

uqGslVectorClass&
uqGslVectorClass::operator*=(const uqGslVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      m_env.fullRank(),
                      "uqGslVectorClass::operator*=()",
                      "different sizes of this and rhs");

  for (unsigned int i = 0; i < size1; ++i) {
    (*this)[i] *= rhs[i];
  }

  return *this;
}

uqGslVectorClass&
uqGslVectorClass::operator/=(const uqGslVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      m_env.fullRank(),
                      "uqGslVectorClass::operator/=()",
                      "different sizes of this and rhs");

  for (unsigned int i = 0; i < size1; ++i) {
    (*this)[i] /= rhs[i];
  }

  return *this;
}

uqGslVectorClass&
uqGslVectorClass::operator+=(const uqGslVectorClass& rhs)
{
  int iRC;
  iRC = gsl_vector_add(m_vec,rhs.m_vec);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslVectorClass::operator+=()",
                    "failed");
  return *this;
}

uqGslVectorClass&
uqGslVectorClass::operator-=(const uqGslVectorClass& rhs)
{
  int iRC;
  iRC = gsl_vector_sub(m_vec,rhs.m_vec);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslVectorClass::operator-=()",
                    "failed");

  return *this;
}

double&
uqGslVectorClass::operator[](unsigned int i)
{
  return *gsl_vector_ptr(m_vec,i);
}

const double&
uqGslVectorClass::operator[](unsigned int i) const
{
  return *gsl_vector_const_ptr(m_vec,i);
}

void
uqGslVectorClass::copy(const uqGslVectorClass& src)
{
  int iRC;
  iRC = gsl_vector_memcpy(this->m_vec, src.m_vec);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.fullRank(),
                    "uqGslVectorClass::copy()",
                    "failed");

  return;
}

unsigned int
uqGslVectorClass::sizeLocal() const
{
  return m_vec->size;
}

unsigned int
uqGslVectorClass::sizeGlobal() const
{
  return m_vec->size;
}

double
uqGslVectorClass::norm2Sq() const
{
  return scalarProduct(*this,*this);
}

double
uqGslVectorClass::norm2() const
{
  return std::sqrt(this->norm2Sq());
}

double
uqGslVectorClass::sumOfComponents() const
{
  double result = 0.;
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += (*this)[i];
  }

  return result;
}

void
uqGslVectorClass::cwSet(double value)
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = value;
  }

  return;
}

void
uqGslVectorClass::cwSetGaussian(const gsl_rng* rng, double mean, double stdDev)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = mean + gsl_ran_gaussian(rng,stdDev);
  }

  return;
}

void
uqGslVectorClass::cwSetGaussian(const gsl_rng* rng, const uqGslVectorClass& meanVec, const uqGslVectorClass& stdDevVec)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = meanVec[i] + gsl_ran_gaussian(rng,stdDevVec[i]);
  }
  return;
}

void
uqGslVectorClass::cwSetUniform(const gsl_rng* rng, const uqGslVectorClass& aVec, const uqGslVectorClass& bVec)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = aVec[i] + (bVec[i]-aVec[i])*gsl_rng_uniform(rng);
  }
  return;
}

void
uqGslVectorClass::cwInvert()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = 1./(*this)[i];
  }

  return;
}

void
uqGslVectorClass::sort()
{
  gsl_sort_vector(m_vec);

  return;
}

void
uqGslVectorClass::print(std::ostream& os) const
{
  unsigned int size = this->sizeLocal();

  if (m_printHorizontally) {
    for (unsigned int i = 0; i < size; ++i) {
      os << (*this)[i]
         << " ";
    }
  }
  else {
    for (unsigned int i = 0; i < size; ++i) {
      os << (*this)[i]
         << std::endl;
    }
  }

  return;
}

gsl_vector*
uqGslVectorClass::data() const
{
  return m_vec;
}

bool
uqGslVectorClass::atLeastOneComponentSmallerThan(const uqGslVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.fullRank(),
                      "uqGslVectorClass::atLeastOneComponentSmallerThan()",
                      "vectors have different sizes");

  bool result = false;
  unsigned int i = 0;
  unsigned int size = this->sizeLocal();
  while ((i < size) && (result == false)) {
    result = ( (*this)[i] < rhs[i] );
    i++;
  };

  return result;
}

bool
uqGslVectorClass::atLeastOneComponentBiggerThan (const uqGslVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.fullRank(),
                      "uqGslVectorClass::atLeastOneComponentBiggerThan()",
                      "vectors have different sizes");

  bool result = false;
  unsigned int i = 0;
  unsigned int size = this->sizeLocal();
  while ((i < size) && (result == false)) {
    result = ( (*this)[i] > rhs[i] );
    i++;
  };

  return result;
}

std::ostream&
operator<<(std::ostream& os, const uqGslVectorClass& obj)
{
  obj.print(os);

  return os;
}

uqGslVectorClass operator/(double a, const uqGslVectorClass& x)
{
  uqGslVectorClass answer(x);
  answer.cwInvert();
  answer *= a;

  return answer;
}

uqGslVectorClass operator/(const uqGslVectorClass& x, const uqGslVectorClass& y)
{
  uqGslVectorClass answer(x);
  answer /= y;

  return answer;
}

uqGslVectorClass operator*(double a, const uqGslVectorClass& x)
{
  uqGslVectorClass answer(x);
  answer *= a;

  return answer;
}

uqGslVectorClass operator*(const uqGslVectorClass& x, const uqGslVectorClass& y)
{
  uqGslVectorClass answer(x);
  answer *= y;

  return answer;
}

double scalarProduct(const uqGslVectorClass& x, const uqGslVectorClass& y)
{
  unsigned int size1 = x.sizeLocal();
  unsigned int size2 = y.sizeLocal();
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      x.env().fullRank(),
                      "scalarProduct()",
                      "different sizes of x and y");

  double result = 0.;
  for (unsigned int i = 0; i < size1; ++i) {
    result += x[i]*y[i];
  }

  return result;
}

uqGslVectorClass operator+(const uqGslVectorClass& x, const uqGslVectorClass& y)
{
  uqGslVectorClass answer(x);
  answer += y;

  return answer;
}

uqGslVectorClass operator-(const uqGslVectorClass& x, const uqGslVectorClass& y)
{
  uqGslVectorClass answer(x);
  answer -= y;

  return answer;
}

