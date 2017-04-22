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

#include <queso/Defines.h>
#include <queso/GslVector.h>
#include <queso/RngBase.h>
#include <algorithm>
#include <gsl/gsl_sort_vector.h>
#include <cmath>

namespace QUESO {

GslVector::GslVector(const BaseEnvironment& env, const Map& map)
  :
  Vector(env,map),
  m_vec        (gsl_vector_calloc(map.NumGlobalElements()))
{
  //std::cout << "Entering GslVector::constructor(1)" << std::endl;

  queso_require_msg(m_vec, "null vector generated");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) map.NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) map.NumGlobalElements(), "incompatible global vec size");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) m_map.NumMyElements(), "incompatible own vec size");

  //std::cout << "In GslVector::constructor(env,map)"
  //          << "\n  m_vec->size             = " << m_vec->size
  //          << "\n  map.NumGlobalElements() = " << map.NumGlobalElements()
  //          << "\n  map.NumMyElements()     = " << map.NumMyElements()
  //          << std::endl;

}

GslVector::GslVector(const BaseEnvironment& env, const Map& map, double value)
  :
  Vector(env,map),
  m_vec        (gsl_vector_calloc(map.NumGlobalElements()))
{
  //std::cout << "Entering GslVector::constructor(2)" << std::endl;

  queso_require_msg(m_vec, "null vector generated");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) map.NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) map.NumGlobalElements(), "incompatible global vec size");

  this->cwSet(value);

  queso_require_equal_to_msg(m_vec->size, (unsigned int) m_map.NumMyElements(), "incompatible own vec size");

  //std::cout << "Leaving GslVector::constructor(2)" << std::endl;
}

GslVector::GslVector(const BaseEnvironment& env, double d1, double d2, const Map& map)
  :
  Vector(env,map),
  m_vec        (gsl_vector_calloc(map.NumGlobalElements()))
{
  //std::cout << "Entering GslVector::constructor(3)" << std::endl;

  queso_require_msg(m_vec, "null vector generated");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) map.NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) map.NumGlobalElements(), "incompatible global vec size");

  for (unsigned int i = 0; i < m_vec->size; ++i) {
    double alpha = (double) i / ((double) m_vec->size - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  queso_require_equal_to_msg(m_vec->size, (unsigned int) m_map.NumMyElements(), "incompatible own vec size");

  //std::cout << "Leaving GslVector::constructor(3)" << std::endl;
}

GslVector::GslVector(const GslVector& v, double start, double end)
  :
  Vector(v.env(),v.map()),
  m_vec        (gsl_vector_calloc(v.sizeLocal()))
{
  //std::cout << "Entering GslVector::constructor(4)" << std::endl;

  queso_require_msg(m_vec, "null vector generated");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) v.map().NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) v.map().NumGlobalElements(), "incompatible global vec size");

  for (unsigned int i = 0; i < m_vec->size; ++i) {
    double alpha = (double) i / ((double) m_vec->size - 1.);
    (*this)[i] = (1. - alpha) * start + alpha * end;
  }

  queso_require_equal_to_msg(m_vec->size, (unsigned int) m_map.NumMyElements(), "incompatible own vec size");

  //std::cout << "Leaving GslVector::constructor(4)" << std::endl;
}

GslVector::GslVector(const GslVector& v)  // mox
  :
  Vector(v.env(),v.map()),
  m_vec        (gsl_vector_calloc(v.sizeLocal()))
{
  //std::cout << "Entering GslVector::constructor(5)" << std::endl;

  // prudenci 2010-06-17 mox
  queso_require_msg(m_vec, "null vector generated");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) v.map().NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec->size, (unsigned int) v.map().NumGlobalElements(), "incompatible global vec size");

  this->copy(v);

  queso_require_equal_to_msg(m_vec->size, (unsigned int) m_map.NumMyElements(), "incompatible own vec size");

  //std::cout << "Leaving GslVector::constructor(5)" << std::endl;
}

GslVector::~GslVector()
{
  if (m_vec) gsl_vector_free(m_vec);
}

GslVector&
GslVector::operator=(const GslVector& rhs)
{
  //std::cout << "In GslVector::operator=()" // mox
  //          << ": setting size1"
  //          << std::endl;
  unsigned int size1 = this->sizeLocal();
  //std::cout << "In GslVector::operator=()" // mox
  //          << ": setting size2"
  //          << std::endl;
  unsigned int size2 = rhs.sizeLocal();
  queso_require_equal_to_msg(size1, size2, "sizes are not compatible");
  this->copy(rhs);
  return *this;
}

GslVector&
GslVector::operator*=(double a)
{
  int iRC;
  iRC = gsl_vector_scale(m_vec,a);
  queso_require_msg(!(iRC), "failed");
  return *this;
}

GslVector&
GslVector::operator/=(double a)
{
  *this *= (1./a);

  return *this;
}

GslVector&
GslVector::operator*=(const GslVector& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  queso_require_equal_to_msg(size1, size2, "different sizes of this and rhs");

  for (unsigned int i = 0; i < size1; ++i) {
    (*this)[i] *= rhs[i];
  }

  return *this;
}

GslVector&
GslVector::operator/=(const GslVector& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  queso_require_equal_to_msg(size1, size2, "different sizes of this and rhs");

  for (unsigned int i = 0; i < size1; ++i) {
    (*this)[i] /= rhs[i];
  }

  return *this;
}

GslVector&
GslVector::operator+=(const GslVector& rhs)
{
  int iRC;
  iRC = gsl_vector_add(m_vec,rhs.m_vec);
  queso_require_msg(!(iRC), "failed");
  return *this;
}

GslVector&
GslVector::operator-=(const GslVector& rhs)
{
  int iRC;
  iRC = gsl_vector_sub(m_vec,rhs.m_vec);
  queso_require_msg(!(iRC), "failed");

  return *this;
}

void
GslVector::copy(const GslVector& src)
{
  this->Vector::base_copy(src);
  int iRC;
  iRC = gsl_vector_memcpy(this->m_vec, src.m_vec);
  queso_require_msg(!(iRC), "failed");

  return;
}

unsigned int
GslVector::sizeLocal() const
{
  // mox
  //std::cout << "Entering GslVector::sizeLocal()"
  //          << ": &m_map = "                << &m_map
  //          << std::endl;
  //std::cout << ", m_map.NumMyElements() = " << m_map.NumMyElements()
  //          << std::endl;
  //std::cout << ", m_vec = "                 << m_vec
  //          << std::endl;
  //std::cout << ", m_vec->size = "           << m_vec->size
  //          << std::endl;

  queso_require_equal_to_msg(m_vec->size, (unsigned int) m_map.NumMyElements(), "incompatible vec size");

  //std::cout << "Leaving GslVector::sizeLocal()"
  //          << ": m_vec = " << m_vec
  //          << ", m_vec->size = " << m_vec->size
  //          << ", &m_map = " << &m_map
  //          << ", m_map.NumMyElements() = " << m_map.NumMyElements()
  //          << std::endl;

  return m_vec->size;
}

unsigned int
GslVector::sizeGlobal() const
{
  queso_require_equal_to_msg(m_vec->size, (unsigned int) m_map.NumGlobalElements(), "incompatible vec size");

  return m_vec->size;
}

double
GslVector::norm2Sq() const
{
  return scalarProduct(*this,*this);
}

double
GslVector::norm2() const
{
  return std::sqrt(this->norm2Sq());
}

double
GslVector::norm1() const
{
  double result = 0.;

  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += fabs((*this)[i]);
  }

  return result;
}

double
GslVector::normInf() const
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
GslVector::sumOfComponents() const
{
  double result = 0.;
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += (*this)[i];
  }

  return result;
}

void
GslVector::cwSet(double value)
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = value;
  }

  return;
}

void
GslVector::cwSetGaussian(double mean, double stdDev)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = mean + m_env.rngObject()->gaussianSample(stdDev);
  }

  return;
}

void
GslVector::cwSetGaussian(const GslVector& meanVec, const GslVector& stdDevVec)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = meanVec[i] + m_env.rngObject()->gaussianSample(stdDevVec[i]);
  }
  return;
}

void
GslVector::cwSetUniform(const GslVector& aVec, const GslVector& bVec)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = aVec[i] + (bVec[i]-aVec[i])*m_env.rngObject()->uniformSample();
  }
  return;
}

void
GslVector::cwSetBeta(const GslVector& alpha, const GslVector& beta)
{
  queso_require_equal_to_msg(this->sizeLocal(), alpha.sizeLocal(), "incompatible alpha size");

  queso_require_equal_to_msg(this->sizeLocal(), beta.sizeLocal(), "incompatible beta size");

  double tmpSample = 0.;
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    tmpSample = m_env.rngObject()->betaSample(alpha[i],beta[i]);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In GslVector::cwSetBeta()"
                              << ": fullRank "   << m_env.fullRank()
                              << ", i = "        << i
                              << ", alpha[i] = " << alpha[i]
                              << ", beta[i] = "  << beta[i]
                              << ", sample = "   << tmpSample
                              << std::endl;
    }
    if ((alpha[i] == 1. ) &&
        (beta [i] == 0.1)) {
      if (tmpSample == 1.) {
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
          *m_env.subDisplayFile() << "Hitting 'sampe = 1' in GslVector::cwSetBeta()"
                                  << ": fullRank "   << m_env.fullRank()
                                  << ", i = "        << i
                                  << ", alpha[i] = " << alpha[i]
                                  << ", beta[i] = "  << beta[i]
                                  << ", sample = "   << tmpSample
                                  << std::endl;
        }
#if 1
        std::cerr << "Hitting 'sample = 1' in GslVector::cwSetBeta()"
                  << ": fullRank "   << m_env.fullRank()
                  << ", i = "        << i
                  << ", alpha[i] = " << alpha[i]
                  << ", beta[i] = "  << beta[i]
                  << ", sample = "   << tmpSample
                  << std::endl;
        do {
          tmpSample = m_env.rngObject()->betaSample(alpha[i],beta[i]);
        } while (tmpSample == 1.);
        std::cerr << "Code was able to get 'sample != 1' in GslVector::cwSetBeta()"
                  << ": fullRank "   << m_env.fullRank()
                  << ", i = "        << i
                  << ", alpha[i] = " << alpha[i]
                  << ", beta[i] = "  << beta[i]
                  << ", sample = "   << tmpSample
                  << std::endl;
      }
#endif
    }
    (*this)[i] = tmpSample;
  }
  return;
}

void
GslVector::cwSetGamma(const GslVector& a, const GslVector& b)
{
  queso_require_equal_to_msg(this->sizeLocal(), a.sizeLocal(), "incompatible a size");

  queso_require_equal_to_msg(this->sizeLocal(), b.sizeLocal(), "incompatible b size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = m_env.rngObject()->gammaSample(a[i],b[i]);
  }
  return;
}

void
GslVector::cwSetInverseGamma(const GslVector& alpha, const GslVector& beta)
{
  queso_require_equal_to_msg(this->sizeLocal(), alpha.sizeLocal(), "incompatible alpha size");

  queso_require_equal_to_msg(this->sizeLocal(), beta.sizeLocal(), "incompatible beta size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = 1./m_env.rngObject()->gammaSample(alpha[i],1./beta[i]);
  }
  return;
}

void
GslVector::cwSetConcatenated(const GslVector& v1, const GslVector& v2)
{
  queso_require_equal_to_msg(this->sizeLocal(), v1.sizeLocal() + v2.sizeLocal(), "incompatible vector sizes");

  for (unsigned int i = 0; i < v1.sizeLocal(); ++i) {
    (*this)[i] = v1[i];
  }

  for (unsigned int i = 0; i < v2.sizeLocal(); ++i) {
    (*this)[v1.sizeLocal()+i] = v2[i];
  }

  return;
}

void
GslVector::cwSetConcatenated(const std::vector<const GslVector* >& vecs)
{
  unsigned int cummulativeSize = 0;
  for (unsigned int i = 0; i < vecs.size(); ++i) {
    GslVector tmpVec(*(vecs[i]));
    for (unsigned int j = 0; j < vecs[i]->sizeLocal(); ++j) {
      (*this)[cummulativeSize+j] = tmpVec[j];
    }
    cummulativeSize += vecs[i]->sizeLocal();
  }

  queso_require_equal_to_msg(this->sizeLocal(), cummulativeSize, "incompatible vector sizes");
  return;
}


void
GslVector::cwSet(unsigned int initialPos, const GslVector& vec)
{
  queso_require_less_msg(initialPos, this->sizeLocal(), "invalid initialPos");

  queso_require_less_equal_msg((initialPos +vec.sizeLocal()), this->sizeLocal(), "invalid vec.sizeLocal()");

  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    (*this)[initialPos+i] = vec[i];
  }

  return;
}

void
GslVector::cwExtract(unsigned int initialPos, GslVector& vec) const
{
  queso_require_less_msg(initialPos, this->sizeLocal(), "invalid initialPos");

  queso_require_less_equal_msg((initialPos +vec.sizeLocal()), this->sizeLocal(), "invalid vec.sizeLocal()");

  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    vec[i] = (*this)[initialPos+i];
  }

  return;
}

void
GslVector::cwInvert()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = 1./(*this)[i];
  }

  return;
}

void
GslVector::cwSqrt()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = sqrt((*this)[i]);
  }

  return;
}

void
GslVector::matlabDiff(
  unsigned int      firstPositionToStoreDiff,
  double            valueForRemainderPosition,
  GslVector& outputVec) const
{
  unsigned int size = this->sizeLocal();

  queso_require_less_equal_msg(firstPositionToStoreDiff, 1, "invalid firstPositionToStoreDiff");

  queso_require_equal_to_msg(size, outputVec.sizeLocal(), "invalid size of outputVecs");

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
GslVector::matlabLinearInterpExtrap(
  const GslVector& x1Vec,
  const GslVector& y1Vec,
  const GslVector& x2Vec)
{
  queso_require_greater_msg(x1Vec.sizeLocal(), 1, "invalid 'x1' size");

  queso_require_equal_to_msg(x1Vec.sizeLocal(), y1Vec.sizeLocal(), "invalid 'x1' and 'y1' sizes");

  queso_require_equal_to_msg(x2Vec.sizeLocal(), this->sizeLocal(), "invalid 'x2' and 'this' sizes");

  for (unsigned int i = 1; i < x1Vec.sizeLocal(); ++i) { // Yes, '1'
    queso_require_greater_msg(x1Vec[i], x1Vec[i-1], "invalid 'x1' values");
  }

  for (unsigned int id2 = 0; id2 < x2Vec.sizeLocal(); ++id2) {
    double x2 = x2Vec[id2];
    unsigned int id1 = 0;
    bool found1 = false;
    for (id1 = 0; id1 < x1Vec.sizeLocal(); ++id1) {
      if (x2 <= x1Vec[id1]) {
        found1 = true;
        break;
      }
    }
    bool makeLinearModel = false;
    double xa = 0.;
    double xb = 0.;
    double ya = 0.;
    double yb = 0.;
    if (x2 == x1Vec[id1]) {
      (*this)[id2] = y1Vec[id1];
    }
    else if (x2 < x1Vec[0]) {
      // Extrapolation case
      makeLinearModel = true;
      xa = x1Vec[0];
      xb = x1Vec[1];
      ya = y1Vec[0];
      yb = y1Vec[1];
    }
    else if (found1 == true) {
      // Interpolation case
      makeLinearModel = true;
      xa = x1Vec[id1-1];
      xb = x1Vec[id1];
      ya = y1Vec[id1-1];
      yb = y1Vec[id1];
    }
    else {
      // Extrapolation case
      makeLinearModel = true;
      xa = x1Vec[x1Vec.sizeLocal()-2];
      xb = x1Vec[x1Vec.sizeLocal()-1];
      ya = y1Vec[x1Vec.sizeLocal()-2];
      yb = y1Vec[x1Vec.sizeLocal()-1];
    }

    if (makeLinearModel) {
      double rate = (yb-ya)/(xb-xa);
      (*this)[id2] = ya + (x2-xa)*rate;
    }
  }



  return;
}

void
GslVector::sort()
{
  gsl_sort_vector(m_vec);

  return;
}

void
GslVector::mpiBcast(int srcRank, const MpiComm& bcastComm)
{
  // Filter out those nodes that should not participate
  if (bcastComm.MyPID() < 0) return;

  // Check 'srcRank'
  queso_require_msg(!((srcRank < 0) || (srcRank >= bcastComm.NumProc())), "invalud srcRank");

  // Check number of participant nodes
  double localNumNodes = 1.;
  double totalNumNodes = 0.;
  bcastComm.Allreduce<double>(&localNumNodes, &totalNumNodes, (int) 1, RawValue_MPI_SUM,
                      "GslVector::mpiBcast()",
                      "failed MPI.Allreduce() for numNodes");
  queso_require_equal_to_msg(((int) totalNumNodes), bcastComm.NumProc(), "inconsistent numNodes");

  // Check that all participant nodes have the same vector size
  double localVectorSize  = this->sizeLocal();
  double sumOfVectorSizes = 0.;
  bcastComm.Allreduce<double>(&localVectorSize, &sumOfVectorSizes, (int) 1, RawValue_MPI_SUM,
                      "GslVector::mpiBcast()",
                      "failed MPI.Allreduce() for vectorSize");

  if ( ((unsigned int) sumOfVectorSizes) != ((unsigned int)(totalNumNodes*localVectorSize)) ) {
    std::cerr << "rank "                 << bcastComm.MyPID()
              << ": sumOfVectorSizes = " << sumOfVectorSizes
              << ", totalNumNodes = "    << totalNumNodes
              << ", localVectorSize = "  << localVectorSize
              << std::endl;
  }
  bcastComm.Barrier();
  queso_require_equal_to_msg(((unsigned int) sumOfVectorSizes), ((unsigned int)(totalNumNodes*localVectorSize)), "inconsistent vectorSize");

  // Ok, bcast data
  std::vector<double> dataBuffer((unsigned int) localVectorSize, 0.);
  if (bcastComm.MyPID() == srcRank) {
    for (unsigned int i = 0; i < dataBuffer.size(); ++i) {
      dataBuffer[i] = (*this)[i];
    }
  }

  bcastComm.Bcast((void *) &dataBuffer[0], (int) localVectorSize, RawValue_MPI_DOUBLE, srcRank,
                  "GslVector::mpiBcast()",
                  "failed MPI.Bcast()");

  if (bcastComm.MyPID() != srcRank) {
    for (unsigned int i = 0; i < dataBuffer.size(); ++i) {
      (*this)[i] = dataBuffer[i];
    }
  }

  return;
}

void
GslVector::mpiAllReduce(RawType_MPI_Op mpiOperation, const MpiComm& opComm, GslVector& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  unsigned int size = this->sizeLocal();
  queso_require_equal_to_msg(size, resultVec.sizeLocal(), "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double srcValue = (*this)[i];
    double resultValue = 0.;
    opComm.Allreduce<double>(&srcValue, &resultValue, (int) 1, mpiOperation,
                     "GslVector::mpiAllReduce()",
                     "failed MPI.Allreduce()");
    resultVec[i] = resultValue;
  }

  return;
}

void
GslVector::mpiAllQuantile(double probability, const MpiComm& opComm, GslVector& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  queso_require_msg(!((probability < 0.) || (1. < probability)), "invalid input");

  unsigned int size = this->sizeLocal();
  queso_require_equal_to_msg(size, resultVec.sizeLocal(), "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double auxDouble = (int) (*this)[i];
    std::vector<double> vecOfDoubles(opComm.NumProc(),0.);
    opComm.Gather<double>(&auxDouble, 1, &vecOfDoubles[0], (int) 1, 0,
                  "GslVector::mpiAllQuantile()",
                  "failed MPI.Gather()");

    std::sort(vecOfDoubles.begin(), vecOfDoubles.end());

    double result = vecOfDoubles[(unsigned int)( probability*((double)(vecOfDoubles.size()-1)) )];

    opComm.Bcast((void *) &result, (int) 1, RawValue_MPI_DOUBLE, 0,
                 "GslVector::mpiAllQuantile()",
                 "failed MPI.Bcast()");

    resultVec[i] = result;
  }

  return;
}

void
GslVector::print(std::ostream& os) const
{
  //std::cout << "In GslVector::print(): before sizelocal()"
  //          << std::endl;
  unsigned int size = this->sizeLocal();
  //std::cout << "In GslVector::print(): after sizelocal()"
  //          << std::endl;

  //std::cout << "In GslVector::print(): before os.flags()"
  //          << std::endl;
  std::ostream::fmtflags curr_fmt = os.flags();
  //std::cout << "In GslVector::print(): after os.flags()"
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
      //std::cout << "In GslVector::print(): where expected"
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

  //std::cout << "In GslVector::print(): before os.flags(curr_fmt)"
  //          << std::endl;
  os.flags(curr_fmt);
  //std::cout << "In GslVector::print(): after os.flags(curr_fmt)"
  //          << std::endl;

  return;
}

void
GslVector::subWriteContents(
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
GslVector::subReadContents(
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
    queso_require_equal_to_msg(tmpString, std::string("="), std::string("string should be the '=' sign"));

    // Read 'zeros(n_positions,n_params)' string
    *filePtrSet.ifsVar >> tmpString;
    //std::cout << "Just read '" << tmpString << "'" << std::endl;
    unsigned int posInTmpString = 6;

    // Isolate 'n_positions' in a string
    char nPositionsString[tmpString.size()-posInTmpString+1];
    unsigned int posInPositionsString = 0;
    do {
      queso_require_less_msg(posInTmpString, tmpString.size(), "symbol ',' not found in first line of file");
      nPositionsString[posInPositionsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ',');
    nPositionsString[posInPositionsString] = '\0';

    // Isolate 'n_params' in a string
    posInTmpString++; // Avoid reading ',' char
    char nParamsString[tmpString.size()-posInTmpString+1];
    unsigned int posInParamsString = 0;
    do {
      queso_require_less_msg(posInTmpString, tmpString.size(), "symbol ')' not found in first line of file");
      nParamsString[posInParamsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ')');
    nParamsString[posInParamsString] = '\0';

    // Convert 'n_positions' and 'n_params' strings to numbers
    unsigned int sizeOfVecInFile = (unsigned int) strtod(nPositionsString,NULL);
    unsigned int numParamsInFile = (unsigned int) strtod(nParamsString,   NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GslVector::subReadContents()"
                              << ": fullRank "            << m_env.fullRank()
                              << ", sizeOfVecInFile = "   << sizeOfVecInFile
                              << ", numParamsInFile = "   << numParamsInFile
                              << ", this->sizeLocal() = " << this->sizeLocal()
                              << std::endl;
    }

    // Check if [size of vec in file] >= [requested sub vec size]
    queso_require_greater_equal_msg(sizeOfVecInFile, subReadSize, "size of vec in file is not big enough");

    // Check if [num params in file] == [num params in current vec]
    queso_require_equal_to_msg(numParamsInFile, numParams, "number of parameters of vec in file is different than number of parameters in this vec object");

    // Code common to any core in a communicator
    unsigned int maxCharsPerLine = 64*numParams; // Up to about 60 characters to represent each parameter value

    unsigned int lineId = 0;
    while (lineId < idOfMyFirstLine) {
      filePtrSet.ifsVar->ignore(maxCharsPerLine,'\n');
      lineId++;
    };

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In GslVector::subReadContents()"
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
      *m_env.subDisplayFile() << "In GslVector::subReadContents()"
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
GslVector::data() const
{
  return m_vec;
}

bool
GslVector::atLeastOneComponentSmallerThan(const GslVector& rhs) const
{
  queso_require_equal_to_msg(this->sizeLocal(), rhs.sizeLocal(), "vectors have different sizes");

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
GslVector::atLeastOneComponentBiggerThan(const GslVector& rhs) const
{
  queso_require_equal_to_msg(this->sizeLocal(), rhs.sizeLocal(), "vectors have different sizes");

  bool result = false;
  unsigned int i = 0;
  unsigned int size = this->sizeLocal();
  while ((i < size) && (result == false)) {
    result = ( (*this)[i] > rhs[i] );
    i++;
  };

  return result;
}

bool
GslVector::atLeastOneComponentSmallerOrEqualThan(const GslVector& rhs) const
{
  queso_require_equal_to_msg(this->sizeLocal(), rhs.sizeLocal(), "vectors have different sizes");

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
GslVector::atLeastOneComponentBiggerOrEqualThan(const GslVector& rhs) const
{
  queso_require_equal_to_msg(this->sizeLocal(), rhs.sizeLocal(), "vectors have different sizes");

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
GslVector::getMaxValue( ) const
{
  return gsl_vector_max( m_vec );
}

double
GslVector::getMinValue( ) const
{
  return gsl_vector_min( m_vec );
}

int
GslVector::getMaxValueIndex( ) const
{
  return gsl_vector_max_index( m_vec );
}

int
GslVector::getMinValueIndex( ) const
{
  return gsl_vector_min_index( m_vec );
}

void
GslVector::getMaxValueAndIndex( double& max_value, int& max_value_index )
{
  max_value = this->getMaxValue();
  max_value_index = this->getMaxValueIndex();

  return;
}

void
GslVector::getMinValueAndIndex( double& min_value, int& min_value_index )
{
  min_value = this->getMinValue();
  min_value_index = this->getMinValueIndex();

  return;
}

GslVector
GslVector::abs() const
{
  GslVector abs_of_this_vec( *this );

  unsigned int size = abs_of_this_vec.sizeLocal();

  for( unsigned int i = 0; i < size; ++i )
    {
      abs_of_this_vec[i] = std::fabs( (*this)[i] );
    }

  return abs_of_this_vec;
}

std::ostream&
operator<<(std::ostream& os, const GslVector& obj)
{
  obj.print(os);

  return os;
}

GslVector operator/(double a, const GslVector& x)
{
  GslVector answer(x);
  answer.cwInvert();
  answer *= a;

  return answer;
}

GslVector operator/(const GslVector& x, const GslVector& y)
{
  GslVector answer(x);
  answer /= y;

  return answer;
}

GslVector operator*(double a, const GslVector& x)
{
  GslVector answer(x);
  answer *= a;

  return answer;
}

GslVector operator*(const GslVector& x, const GslVector& y)
{
  GslVector answer(x);
  answer *= y;

  return answer;
}

double scalarProduct(const GslVector& x, const GslVector& y)
{
  unsigned int size1 = x.sizeLocal();
  unsigned int size2 = y.sizeLocal();
  queso_require_equal_to_msg(size1, size2, "different sizes of x and y");

  double result = 0.;
  for (unsigned int i = 0; i < size1; ++i) {
    result += x[i]*y[i];
  }

  return result;
}

GslVector operator+(const GslVector& x, const GslVector& y)
{
  GslVector answer(x);
  answer += y;

  return answer;
}

GslVector operator-(const GslVector& x, const GslVector& y)
{
  GslVector answer(x);
  answer -= y;

  return answer;
}

bool
operator== (const GslVector& lhs, const GslVector& rhs)
{
  bool answer = true;

  unsigned int size1 = lhs.sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  queso_require_equal_to_msg(size1, size2, "different sizes of lhs and rhs");

  for (unsigned int i = 0; i < size1; ++i) {
    if (lhs[i] != rhs[i]) {
      answer = false;
      break;
    }
  }

  return answer;
}

}  // End namespace QUESO
