//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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

#include <uqGslVector.h>
#include <uqDefines.h>
#include <gsl/gsl_sort_vector.h>
#include <cmath>

uqGslVectorClass::uqGslVectorClass()
  :
  uqVectorClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(), default",
                      "should not be used by user");
}

uqGslVectorClass::uqGslVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map)
  :
  uqVectorClass(env,map),
  m_vec        (gsl_vector_calloc(map.NumGlobalElements()))
{
  //std::cout << "Entering uqGslVectorClass::constructor(1)" << std::endl;

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(1)",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(1)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) map.NumGlobalElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(1)",
                      "incompatible global vec size");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(1)",
                      "incompatible own vec size");

  //std::cout << "In uqGslVectorClass::constructor(env,map)"
  //          << "\n  m_vec->size             = " << m_vec->size
  //          << "\n  map.NumGlobalElements() = " << map.NumGlobalElements()
  //          << "\n  map.NumMyElements()     = " << map.NumMyElements()
  //          << std::endl;

  //std::cout << "Leaving uqGslVectorClass::constructor(1)" << std::endl;
}

uqGslVectorClass::uqGslVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double value)
  :
  uqVectorClass(env,map),
  m_vec        (gsl_vector_calloc(map.NumGlobalElements()))
{
  //std::cout << "Entering uqGslVectorClass::constructor(2)" << std::endl;

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(2)",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(2)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) map.NumGlobalElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(2)",
                      "incompatible global vec size");

  this->cwSet(value);

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(2)",
                      "incompatible own vec size");

  //std::cout << "Leaving uqGslVectorClass::constructor(2)" << std::endl;
}

uqGslVectorClass::uqGslVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const uqMapClass& map)
  :
  uqVectorClass(env,map),
  m_vec        (gsl_vector_calloc(map.NumGlobalElements()))
{
  //std::cout << "Entering uqGslVectorClass::constructor(3)" << std::endl;

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(3), linspace",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(3)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) map.NumGlobalElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(3)",
                      "incompatible global vec size");

  for (unsigned int i = 0; i < m_vec->size; ++i) {
    double alpha = (double) i / ((double) m_vec->size - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(3)",
                      "incompatible own vec size");

  //std::cout << "Leaving uqGslVectorClass::constructor(3)" << std::endl;
}

uqGslVectorClass::uqGslVectorClass(const uqGslVectorClass& v, double d1, double d2)
  :
  uqVectorClass(v.env(),v.map()),
  m_vec        (gsl_vector_calloc(v.sizeLocal()))
{
  //std::cout << "Entering uqGslVectorClass::constructor(4)" << std::endl;

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(4), linspace",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) v.map().NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(4)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) v.map().NumGlobalElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(4)",
                      "incompatible global vec size");

  for (unsigned int i = 0; i < m_vec->size; ++i) {
    double alpha = (double) i / ((double) m_vec->size - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(4)",
                      "incompatible own vec size");

  //std::cout << "Leaving uqGslVectorClass::constructor(4)" << std::endl;
}

uqGslVectorClass::uqGslVectorClass(const uqGslVectorClass& v)  // mox
  :
  uqVectorClass(v.env(),v.map()),
  m_vec        (gsl_vector_calloc(v.sizeLocal()))
{
  //std::cout << "Entering uqGslVectorClass::constructor(5)" << std::endl;

  // prudenci 2010-06-17 mox
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(5), copy",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) v.map().NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(5)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) v.map().NumGlobalElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(5)",
                      "incompatible global vec size");

  this->copy(v);

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::constructor(5)",
                      "incompatible own vec size");

  //std::cout << "Leaving uqGslVectorClass::constructor(5)" << std::endl;
}

uqGslVectorClass::~uqGslVectorClass()
{
  if (m_vec) gsl_vector_free(m_vec);
}

uqGslVectorClass&
uqGslVectorClass::operator=(const uqGslVectorClass& rhs)
{
  //std::cout << "In uqGslVectorClass::operator=()" // mox
  //          << ": setting size1"
  //          << std::endl;
  unsigned int size1 = this->sizeLocal();
  //std::cout << "In uqGslVectorClass::operator=()" // mox
  //          << ": setting size2"
  //          << std::endl;
  unsigned int size2 = rhs.sizeLocal();
  UQ_FATAL_TEST_MACRO(size1 != size2, // mox
                      m_env.worldRank(),
                      "uqGslVectorClass::operator=()",
                      "sizes are not compatible");
  this->copy(rhs);
  return *this;
}

uqGslVectorClass&
uqGslVectorClass::operator*=(double a)
{
  int iRC;
  iRC = gsl_vector_scale(m_vec,a);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
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
                      m_env.worldRank(),
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
                      m_env.worldRank(),
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
                    m_env.worldRank(),
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
                    m_env.worldRank(),
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
  this->uqVectorClass::copy(src); // prudenci 2010-06-17 mox
  int iRC;
  iRC = gsl_vector_memcpy(this->m_vec, src.m_vec);
  UQ_FATAL_RC_MACRO(iRC,
                    m_env.worldRank(),
                    "uqGslVectorClass::copy()",
                    "failed");

  return;
}

unsigned int
uqGslVectorClass::sizeLocal() const
{
  // mox
  //std::cout << "Entering uqGslVectorClass::sizeLocal()"
  //          << ": &m_map = "                << &m_map
  //          << std::endl;
  //std::cout << ", m_map.NumMyElements() = " << m_map.NumMyElements()
  //          << std::endl;
  //std::cout << ", m_vec = "                 << m_vec
  //          << std::endl;
  //std::cout << ", m_vec->size = "           << m_vec->size
  //          << std::endl;

  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::sizeLocal()",
                      "incompatible vec size");

  //std::cout << "Leaving uqGslVectorClass::sizeLocal()"
  //          << ": m_vec = " << m_vec
  //          << ", m_vec->size = " << m_vec->size
  //          << ", &m_map = " << &m_map
  //          << ", m_map.NumMyElements() = " << m_map.NumMyElements()
  //          << std::endl;

  return m_vec->size;
}

unsigned int
uqGslVectorClass::sizeGlobal() const
{
  UQ_FATAL_TEST_MACRO(m_vec->size != (unsigned int) m_map.NumGlobalElements(),
                      m_env.worldRank(),
                      "uqGslVectorClass::sizeGlobal()",
                      "incompatible vec size");

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
uqGslVectorClass::norm1() const
{
  double result = 0.;

  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += fabs((*this)[i]);
  }

  return result;
}

double
uqGslVectorClass::normInf() const
{
  double result = 0.;

  unsigned int size = this->sizeLocal();
  double aux = 0.;
  for (unsigned int i = 0; i < size; ++i) {
    aux = fabs((*this)[i]);
    if (aux > result) result = aux;
  }

  return result;
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
uqGslVectorClass::cwSetBeta(const gsl_rng* rng, const uqGslVectorClass& alpha, const uqGslVectorClass& beta)
{
  UQ_FATAL_TEST_MACRO(this->sizeLocal() != alpha.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSetBeta()",
                      "incompatible alpha size");

  UQ_FATAL_TEST_MACRO(this->sizeLocal() != beta.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSetBeta()",
                      "incompatible beta size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = gsl_ran_beta(rng,alpha[i],beta[i]);
  }
  return;
}

void
uqGslVectorClass::cwSetGamma(const gsl_rng* rng, const uqGslVectorClass& a, const uqGslVectorClass& b)
{
  UQ_FATAL_TEST_MACRO(this->sizeLocal() != a.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSetGamma()",
                      "incompatible a size");

  UQ_FATAL_TEST_MACRO(this->sizeLocal() != b.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSetGamma()",
                      "incompatible b size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = gsl_ran_gamma(rng,a[i],b[i]);
  }
  return;
}

void
uqGslVectorClass::cwSetInverseGamma(const gsl_rng* rng, const uqGslVectorClass& alpha, const uqGslVectorClass& beta)
{
  UQ_FATAL_TEST_MACRO(this->sizeLocal() != alpha.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSetInverseGamma()",
                      "incompatible alpha size");

  UQ_FATAL_TEST_MACRO(this->sizeLocal() != beta.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSetInverseGamma()",
                      "incompatible beta size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = 1./gsl_ran_gamma(rng,alpha[i],1./beta[i]);
  }
  return;
}

void
uqGslVectorClass::cwSetConcatenated(const uqGslVectorClass& v1, const uqGslVectorClass& v2)
{
  UQ_FATAL_TEST_MACRO(this->sizeLocal() != v1.sizeLocal() + v2.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSetConcatenated(1)",
                      "incompatible vector sizes");

  for (unsigned int i = 0; i < v1.sizeLocal(); ++i) {
    (*this)[i] = v1[i];
  }

  for (unsigned int i = 0; i < v2.sizeLocal(); ++i) {
    (*this)[v1.sizeLocal()+i] = v2[i];
  }

  return;
}

void
uqGslVectorClass::cwSetConcatenated(const std::vector<const uqGslVectorClass*>& vecs)
{
  unsigned int cummulativeSize = 0;
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    uqGslVectorClass tmpVec(*(vecs[i]));
    for (unsigned int j = 0; j < vecs[i]->sizeLocal(); ++j) {
      (*this)[cummulativeSize+j] = tmpVec[j];
    }
    cummulativeSize += vecs[i]->sizeLocal();
  }

  UQ_FATAL_TEST_MACRO(this->sizeLocal() != cummulativeSize,
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSetConcatenated(1)",
                      "incompatible vector sizes");
  return;
}


void
uqGslVectorClass::cwSet(unsigned int initialPos, const uqGslVectorClass& vec)
{
  UQ_FATAL_TEST_MACRO(initialPos >= this->sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSet()",
                      "invalid initialPos");

  UQ_FATAL_TEST_MACRO((initialPos +vec.sizeLocal()) > this->sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwSet()",
                      "invalid vec.sizeLocal()");

  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    (*this)[initialPos+i] = vec[i];
  }

  return;
}

void
uqGslVectorClass::cwExtract(unsigned int initialPos, uqGslVectorClass& vec) const
{
  UQ_FATAL_TEST_MACRO(initialPos >= this->sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwExtract()",
                      "invalid initialPos");

  UQ_FATAL_TEST_MACRO((initialPos +vec.sizeLocal()) > this->sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::cwExtract()",
                      "invalid vec.sizeLocal()");

  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    vec[i] = (*this)[initialPos+i];
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
uqGslVectorClass::cwSqrt()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = sqrt((*this)[i]);
  }

  return;
}

void
uqGslVectorClass::matlabDiff(
  unsigned int      firstPositionToStoreDiff,
  double            valueForRemainderPosition,
  uqGslVectorClass& outputVec) const
{
  unsigned int size = this->sizeLocal();

  UQ_FATAL_TEST_MACRO(firstPositionToStoreDiff > 1,
                      m_env.worldRank(),
                      "uqGslVectorClass::matlabDiff()",
                      "invalid firstPositionToStoreDiff");

  UQ_FATAL_TEST_MACRO(size != outputVec.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::matlabDiff()",
                      "invalid size of outputVecs");

  for (unsigned int i = 0; i < (size-1); ++i) {
    outputVec[firstPositionToStoreDiff+i] = (*this)[i+1]-(*this)[i];
  }
  if (firstPositionToStoreDiff == 0) {
    outputVec[size-1] = valueForRemainderPosition;
  }
  else {
    outputVec[0] = valueForRemainderPosition;
  }

  return;
}

void
uqGslVectorClass::matlabLinearInterpExtrap(
  const uqGslVectorClass& x1Vec,
  const uqGslVectorClass& y1Vec,
  const uqGslVectorClass& x2Vec)
{
  // todo_r

  return;
}

void
uqGslVectorClass::sort()
{
  gsl_sort_vector(m_vec);

  return;
}

void
uqGslVectorClass::mpiBcast(int srcRank, const uqMpiCommClass& bcastComm)
{
  // Filter out those nodes that should not participate
  if (bcastComm.MyPID() < 0) return;

  // Check 'srcRank'
  UQ_FATAL_TEST_MACRO((srcRank < 0) || (srcRank >= bcastComm.NumProc()),
                      m_env.worldRank(),
                      "uqGslVectorClass::mpiBcast()",
                      "invalud srcRank");

  // Check number of participant nodes
  double localNumNodes = 1.;
  double totalNumNodes = 0.;
  bcastComm.Allreduce((void *) &localNumNodes, (void *) &totalNumNodes, (int) 1, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                      "uqGslVectorClass::mpiBcast()",
                      "failed MPI.Allreduce() for numNodes");
  UQ_FATAL_TEST_MACRO(((int) totalNumNodes) != bcastComm.NumProc(),
                      m_env.worldRank(),
                      "uqGslVectorClass::mpiBcast()",
                      "inconsistent numNodes");

  // Check that all participant nodes have the same vector size
  double localVectorSize  = this->sizeLocal();
  double sumOfVectorSizes = 0.; 
  bcastComm.Allreduce((void *) &localVectorSize, (void *) &sumOfVectorSizes, (int) 1, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                      "uqGslVectorClass::mpiBcast()",
                      "failed MPI.Allreduce() for vectorSize");

  if ( ((unsigned int) sumOfVectorSizes) != ((unsigned int)(totalNumNodes*localVectorSize)) ) {
    std::cerr << "rank "                 << bcastComm.MyPID()
              << ": sumOfVectorSizes = " << sumOfVectorSizes
              << ", totalNumNodes = "    << totalNumNodes
              << ", localVectorSize = "  << localVectorSize
              << std::endl;
  }
  bcastComm.Barrier();
  UQ_FATAL_TEST_MACRO(((unsigned int) sumOfVectorSizes) != ((unsigned int)(totalNumNodes*localVectorSize)),
                      m_env.worldRank(),
                      "uqGslVectorClass::mpiBcast()",
                      "inconsistent vectorSize");

  // Ok, bcast data
  std::vector<double> dataBuffer((unsigned int) localVectorSize, 0.);
  if (bcastComm.MyPID() == srcRank) {
    for (unsigned int i = 0; i < dataBuffer.size(); ++i) {
      dataBuffer[i] = (*this)[i];
    }
  }

  bcastComm.Bcast((void *) &dataBuffer[0], (int) localVectorSize, uqRawValue_MPI_DOUBLE, srcRank,
                  "uqGslVectorClass::mpiBcast()",
                  "failed MPI.Bcast()");

  if (bcastComm.MyPID() != srcRank) {
    for (unsigned int i = 0; i < dataBuffer.size(); ++i) {
      (*this)[i] = dataBuffer[i];
    }
  }

  return;
}

void
uqGslVectorClass::mpiAllReduce(uqRawType_MPI_Op mpiOperation, const uqMpiCommClass& opComm, uqGslVectorClass& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  unsigned int size = this->sizeLocal();
  UQ_FATAL_TEST_MACRO(size != resultVec.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::mpiAllReduce()",
                      "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double srcValue = (*this)[i];
    double resultValue = 0.;
    opComm.Allreduce((void *) &srcValue, (void *) &resultValue, (int) 1, uqRawValue_MPI_DOUBLE, mpiOperation,
                     "uqGslVectorClass::mpiAllReduce()",
                     "failed MPI.Allreduce()");
    resultVec[i] = resultValue;
  }

  return;
}

void
uqGslVectorClass::mpiAllQuantile(double probability, const uqMpiCommClass& opComm, uqGslVectorClass& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  UQ_FATAL_TEST_MACRO((probability < 0.) || (1. < probability),
                      m_env.worldRank(),
                      "uqGslVectorClass::mpiAllQuantile()",
                      "invalid input");

  unsigned int size = this->sizeLocal();
  UQ_FATAL_TEST_MACRO(size != resultVec.sizeLocal(),
                      m_env.worldRank(),
                      "uqGslVectorClass::mpiAllQuantile()",
                      "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double auxDouble = (int) (*this)[i];
    std::vector<double> vecOfDoubles(opComm.NumProc(),0.);
    opComm.Gather((void *) &auxDouble, 1, uqRawValue_MPI_DOUBLE, (void *) &vecOfDoubles[0], (int) 1, uqRawValue_MPI_DOUBLE, 0,
                  "uqGslVectorClass::mpiAllQuantile()",
                  "failed MPI.Gather()");

    std::sort(vecOfDoubles.begin(), vecOfDoubles.end());

    double result = vecOfDoubles[(unsigned int)( probability*((double)(vecOfDoubles.size()-1)) )];

    opComm.Bcast((void *) &result, (int) 1, uqRawValue_MPI_DOUBLE, 0,
                 "uqGslVectorClass::mpiAllQuantile()",
                 "failed MPI.Bcast()");

    resultVec[i] = result;
  }

  return;
}

void
uqGslVectorClass::print(std::ostream& os) const
{
  //std::cout << "In uqGslVectorClass::print(): before sizelocal()"
  //          << std::endl;
  unsigned int size = this->sizeLocal();
  //std::cout << "In uqGslVectorClass::print(): after sizelocal()"
  //          << std::endl;

  //std::cout << "In uqGslVectorClass::print(): before os.flags()"
  //          << std::endl;
  std::ostream::fmtflags curr_fmt = os.flags();
  //std::cout << "In uqGslVectorClass::print(): after os.flags()"
  //          << std::endl;

  if (m_printScientific) {
    unsigned int savedPrecision = os.precision();
    os.precision(16);

    if (m_printHorizontally) {
      for (unsigned int i = 0; i < size; ++i) {
        os << std::scientific << (*this)[i]
           << " ";
      }
    }
    else {
      for (unsigned int i = 0; i < size; ++i) {
        os << std::scientific << (*this)[i]
           << std::endl;
      }
    }

    os.precision(savedPrecision);
  }
  else {
    if (m_printHorizontally) {
      //std::cout << "In uqGslVectorClass::print(): where expected"
      //          << std::endl;
      for (unsigned int i = 0; i < size; ++i) {
        os << std::dec << (*this)[i]
           << " ";
      }
    }
    else {
      for (unsigned int i = 0; i < size; ++i) {
        os << std::dec << (*this)[i]
           << std::endl;
      }
    }
  }

  //std::cout << "In uqGslVectorClass::print(): before os.flags(curr_fmt)"
  //          << std::endl;
  os.flags(curr_fmt);
  //std::cout << "In uqGslVectorClass::print(): after os.flags(curr_fmt)"
  //          << std::endl;

  return;
}

void
uqGslVectorClass::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqGslVectorClass::subWriteContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "uqGslVectorClass::subWriteContents()",
                      "implemented just for sequential vectors for now");

  uqFilePtrSetStruct filePtrSet;
  if (m_env.openOutputFile(fileName,
                           fileType, // "m or hdf"
                           allowedSubEnvIds,
                           false,
                           filePtrSet)) {
    *filePtrSet.ofsVar << varNamePrefix << "_sub" << m_env.subIdString() << " = zeros(" << this->sizeLocal()
                       << ","                                                           << 1
                       << ");"
                       << std::endl;
    *filePtrSet.ofsVar << varNamePrefix << "_sub" << m_env.subIdString() << " = [";

    bool savedVectorPrintScientific   = this->getPrintScientific();
    bool savedVectorPrintHorizontally = this->getPrintHorizontally();
    this->setPrintScientific  (true);
    this->setPrintHorizontally(false);
    *filePtrSet.ofsVar << *this;
                     //<< std::endl; // No need for 'endl' because horizontally = 'false'
    this->setPrintHorizontally(savedVectorPrintHorizontally);
    this->setPrintScientific  (savedVectorPrintScientific);

    *filePtrSet.ofsVar << "];\n";

    m_env.closeFile(filePtrSet,fileType);
  }

  return;
}

void
uqGslVectorClass::subReadContents(
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds)
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqGslVectorClass::subReadContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "uqGslVectorClass::subReadContents()",
                      "implemented just for sequential vectors for now");

  uqFilePtrSetStruct filePtrSet;
  if (m_env.openInputFile(fileName,
                          fileType, // "m or hdf"
                          allowedSubEnvIds,
                          filePtrSet)) {
    double subReadSize = this->sizeLocal();

    // In the logic below, the id of a line' begins with value 0 (zero)
    unsigned int idOfMyFirstLine = 1;
    unsigned int idOfMyLastLine = this->sizeLocal();
    unsigned int numParams = 1; // Yes, just '1'

    // Read number of chain positions in the file by taking care of the first line,
    // which resembles something like 'variable_name = zeros(n_positions,m_params);'
    std::string tmpString;

    // Read 'variable name' string
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Just read '" << tmpString << "'" << std::endl;

    // Read '=' sign
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Just read '" << tmpString << "'" << std::endl;
    UQ_FATAL_TEST_MACRO(tmpString != "=",
                        m_env.worldRank(),
                        "uqGslVectorClass::subReadContents()",
                        "string should be the '=' sign");

    // Read 'zeros(n_positions,n_params)' string
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Just read '" << tmpString << "'" << std::endl;
    unsigned int posInTmpString = 6;

    // Isolate 'n_positions' in a string
    char nPositionsString[tmpString.size()-posInTmpString+1];
    unsigned int posInPositionsString = 0;
    do {
      UQ_FATAL_TEST_MACRO(posInTmpString >= tmpString.size(),
                          m_env.worldRank(),
                          "uqGslVectorClass::subReadContents()",
                          "symbol ',' not found in first line of file");
      nPositionsString[posInPositionsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ',');
    nPositionsString[posInPositionsString] = '\0';

    // Isolate 'n_params' in a string
    posInTmpString++; // Avoid reading ',' char
    char nParamsString[tmpString.size()-posInTmpString+1];
    unsigned int posInParamsString = 0;
    do {
      UQ_FATAL_TEST_MACRO(posInTmpString >= tmpString.size(),
                          m_env.worldRank(),
                          "uqGslVectorClass::subReadContents()",
                          "symbol ')' not found in first line of file");
      nParamsString[posInParamsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ')');
    nParamsString[posInParamsString] = '\0';

    // Convert 'n_positions' and 'n_params' strings to numbers
    unsigned int sizeOfVecInFile = (unsigned int) strtod(nPositionsString,NULL);
    unsigned int numParamsInFile = (unsigned int) strtod(nParamsString,   NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGslVectorClass::subReadContents()"
                              << ": fullRank "            << m_env.fullRank()
                              << ", sizeOfVecInFile = "   << sizeOfVecInFile
                              << ", numParamsInFile = "   << numParamsInFile
                              << ", this->sizeLocal() = " << this->sizeLocal()
                              << std::endl;
    }

    // Check if [size of vec in file] >= [requested sub vec size]
    UQ_FATAL_TEST_MACRO(sizeOfVecInFile < subReadSize,
                        m_env.worldRank(),
                        "uqGslVectorClass::subReadContents()",
                        "size of vec in file is not big enough");

    // Check if [num params in file] == [num params in current vec]
    UQ_FATAL_TEST_MACRO(numParamsInFile != numParams,
                        m_env.worldRank(),
                        "uqGslVectorClass::subReadContents()",
                        "number of parameters of vec in file is different than number of parameters in this vec object");

    // Code common to any core in a communicator
    unsigned int maxCharsPerLine = 64*numParams; // Up to about 60 characters to represent each parameter value

    unsigned int lineId = 0;
    while (lineId < idOfMyFirstLine) {
      filePtrSet.ifsVar->ignore(maxCharsPerLine,'\n');
      lineId++;
    };

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGslVectorClass::subReadContents()"
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
                        "uqGslVectorClass::subReadContents()",
                        "in core 0, string should be the '=' sign");

    // Take into account the ' [' portion
    std::streampos tmpPos = filePtrSet.ifsVar->tellg();
    filePtrSet.ifsVar->seekg(tmpPos+(std::streampos)2);

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqGslVectorClass::subReadContents()"
                              << ": beginning to read lines with numbers only"
                              << ", lineId = " << lineId
                              << ", idOfMyFirstLine = " << idOfMyFirstLine
                              << ", idOfMyLastLine = " << idOfMyLastLine
                              << std::endl;
    }

    while (lineId <= idOfMyLastLine) {
      *filePtrSet.ifsVar >> (*this)[lineId - idOfMyFirstLine];
      lineId++;
    };

    m_env.closeFile(filePtrSet,fileType);
  }

  return;
}

gsl_vector*
uqGslVectorClass::data() const
{
  return m_vec;
}

bool
uqGslVectorClass::atLeastOneComponentSmallerOrEqualThan(const uqGslVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslVectorClass::atLeastOneComponentSmallerOrEqualThan()",
                      "vectors have different sizes");

  bool result = false;
  unsigned int i = 0;
  unsigned int size = this->sizeLocal();
  while ((i < size) && (result == false)) {
    result = ( (*this)[i] <= rhs[i] ); // prudencio 2012-02-06
    i++;
  };

  return result;
}

bool
uqGslVectorClass::atLeastOneComponentBiggerOrEqualThan(const uqGslVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "uqGslVectorClass::atLeastOneComponentBiggerOrEqualThan()",
                      "vectors have different sizes");

  bool result = false;
  unsigned int i = 0;
  unsigned int size = this->sizeLocal();
  while ((i < size) && (result == false)) {
    result = ( (*this)[i] >= rhs[i] ); // prudencio 2012-02-06
    i++;
  };

  return result;
}

double
uqGslVectorClass::getMaxValue( ) const
{
  return gsl_vector_max( m_vec );
}

double
uqGslVectorClass::getMinValue( ) const
{
  return gsl_vector_min( m_vec );
}

int
uqGslVectorClass::getMaxValueIndex( ) const
{
  return gsl_vector_max_index( m_vec );
}

int
uqGslVectorClass::getMinValueIndex( ) const
{
  return gsl_vector_min_index( m_vec );
}

void
uqGslVectorClass::getMaxValueAndIndex( double& max_value, int& max_value_index )
{
  max_value = this->getMaxValue();
  max_value_index = this->getMaxValueIndex();

  return;
}

void
uqGslVectorClass::getMinValueAndIndex( double& min_value, int& min_value_index )
{
  min_value = this->getMinValue();
  min_value_index = this->getMinValueIndex();

  return;
}

uqGslVectorClass
uqGslVectorClass::abs() const
{
  uqGslVectorClass abs_of_this_vec( *this );

  unsigned int size = abs_of_this_vec.sizeLocal();

  for( unsigned int i = 0; i < size; ++i )
    {
      abs_of_this_vec[i] = std::fabs( (*this)[i] );
    }

  return abs_of_this_vec;

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
                      x.env().worldRank(),
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

bool
operator== (const uqGslVectorClass& lhs, const uqGslVectorClass& rhs)
{
  bool answer = true;

  unsigned int size1 = lhs.sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      lhs.env().worldRank(),
                      "operator==()",
                      "different sizes of lhs and rhs");

  for (unsigned int i = 0; i < size1; ++i) {
    if (lhs[i] != rhs[i]) {
      answer = false;
      break;
    }
  }

  return answer;
}

