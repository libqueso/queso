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

#ifdef QUESO_HAS_TRILINOS

#include <queso/TeuchosVector.h>
#include <queso/RngBase.h>

namespace QUESO {

using std:: cout;
using std:: endl;

// constructor with dimension -----------------------
TeuchosVector::TeuchosVector(const BaseEnvironment& env, const Map& map)
  :
  Vector(env,map)
{
  m_vec.size(map.NumGlobalElements());

  queso_require_equal_to_msg(m_vec.length(), map.NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec.length(), map.NumGlobalElements(), "incompatible global vec size");

  queso_require_equal_to_msg(m_vec.length(), m_map.NumMyElements(), "incompatible own vec size");

  //std::cout << "In TeuchosVector::constructor(env,map)"
  //          << "\n  m_vec.length()             = " << m_vec.length()
  //          << "\n  map.NumGlobalElements() = " << map.NumGlobalElements()
  //          << "\n  map.NumMyElements()     = " << map.NumMyElements()
  //          << std::endl;
}

// ---------------------------------------------------
TeuchosVector::TeuchosVector(const BaseEnvironment& env, const Map& map, double value)
  :
  Vector(env,map)
{
  m_vec.size(map.NumGlobalElements());
  m_vec = value;

  queso_require_equal_to_msg(m_vec.length(), map.NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec.length(), map.NumGlobalElements(), "incompatible global vec size");

  queso_require_equal_to_msg(m_vec.length(), m_map.NumMyElements(), "incompatible own vec size");
}

// ---------------------------------------------------
TeuchosVector::TeuchosVector(const BaseEnvironment& env, double d1, double d2, const Map& map)
  :
  Vector(env,map)
  {
  m_vec.size(map.NumGlobalElements());

  queso_require_equal_to_msg(m_vec.length(), map.NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec.length(), map.NumGlobalElements(), "incompatible global vec size");

  for (int i = 0; i < m_vec.length(); ++i) {
    double alpha = (double) i / ((double) m_vec.length() - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  queso_require_equal_to_msg(m_vec.length(), m_map.NumMyElements(), "incompatible own vec size");
}

// ---------------------------------------------------
TeuchosVector::TeuchosVector(const TeuchosVector& v, double d1, double d2)
  :
  Vector(v.env(),v.map())
{
  m_vec.size(v.sizeLocal());

  queso_require_equal_to_msg(m_vec.length(), v.map().NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec.length(), v.map().NumGlobalElements(), "incompatible global vec size");

  for ( int i = 0; i < m_vec.length(); ++i) {
    double alpha = (double) i / ((double) m_vec.length() - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  queso_require_equal_to_msg(m_vec.length(), m_map.NumMyElements(), "incompatible own vec size");
}

// ---------------------------------------------------
TeuchosVector::TeuchosVector(const TeuchosVector& v)  // mox
  :
  Vector(v.env(),v.map())
 {
   m_vec.size(v.sizeLocal());

  queso_require_equal_to_msg(m_vec.length(), v.map().NumMyElements(), "incompatible local vec size");

  queso_require_equal_to_msg(m_vec.length(), v.map().NumGlobalElements(), "incompatible global vec size");
  this->copy(v);

  queso_require_equal_to_msg(m_vec.length(), m_map.NumMyElements(), "incompatible own vec size");
}


// destructor --------------------------------------
TeuchosVector::~TeuchosVector()
{
};

// Set methods ------------------------------------
//-------------------------------------------------
// assigns values a to all entrances in the vector
TeuchosVector& TeuchosVector::operator=(double a)
{
  m_vec.putScalar(a);
  return *this;
}

//-------------------------------------------------
TeuchosVector& TeuchosVector::operator=(const TeuchosVector& rhs)
{
  unsigned int size1 = m_vec.length();
  unsigned int size2 = rhs.sizeLocal();

  queso_require_equal_to_msg(size1, size2, "the vectors do NOT have the same size.\n");

  if (size1==size2){
    for (unsigned int i=0;i<size1;i++){
    m_vec[i]=rhs[i];
    }
  }

  return *this;
}

//-------------------------------------------------
TeuchosVector& TeuchosVector::operator*=(double a)
{
  m_vec.scale(a); 	//Scale this vector by a; *this = a*this.
  return *this;
}

 //-------------------------------------------------
TeuchosVector& TeuchosVector::operator/=(double a)
{
  (*this) *= (1.0/a);
  return *this;
}

//-------------------------------------------------
TeuchosVector& TeuchosVector::operator*=(const TeuchosVector& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  queso_require_equal_to_msg(size1, size2, "the vectors do NOT have the same size.\n");
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] *= rhs[i];
    }
 }
 return *this;
}

//-------------------------------------------------
TeuchosVector& TeuchosVector::operator/=(const TeuchosVector& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();

  queso_require_equal_to_msg(size1, size2, "the vectors do NOT have the same size.\n");
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] /= rhs[i];
    }
  }
  return *this;
}

//-------------------------------------------------
TeuchosVector& TeuchosVector::operator+=(const TeuchosVector& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();

  queso_require_equal_to_msg(size1, size2, "the vectors do NOT have the same size.\n");

  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] += rhs[i];
    }
  }
  return *this;
}

//-------------------------------------------------
TeuchosVector& TeuchosVector::operator-=(const TeuchosVector& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();

  queso_require_equal_to_msg(size1, size2, "the vectors do NOT have the same size.\n");

  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] -= rhs[i];
    }
  }

  return *this;
}


// Accessor methods --------------------------------
//-------------------------------------------------
double& TeuchosVector::operator[](unsigned int i)
{
  return m_vec[i];
}

//-------------------------------------------------
const double& TeuchosVector::operator[](unsigned int i) const
{
  return m_vec[i];
}

// Attribute methods ------------------------------
//-------------------------------------------------
unsigned int TeuchosVector::sizeLocal() const
{
  queso_require_equal_to_msg(m_vec.length(), (int) m_map.NumMyElements(), "incompatible vec size");

  return m_vec.length();
}

//-------------------------------------------------
unsigned int TeuchosVector::sizeGlobal() const
{
   queso_require_equal_to_msg(m_vec.length(), (int) m_map.NumGlobalElements(), "incompatible vec size");
  return m_vec.length();
}

//-------------------------------------------------
// TODO: needs to be checked. It may not be used at all. Kemelli 4/30/13.
double*
TeuchosVector::values()
{
  return  m_vec.values();
};

// Getting Max and Min values -------------------
// Max ------------------------------------------
double TeuchosVector::getMaxValue( ) const  //dummy
{
  const unsigned int size = this->sizeLocal();
  std::vector<double> aux;

  for (unsigned int i=0; i<size; i++ ) {
    aux.push_back((*this)[i]) ;
  }

  return *max_element (aux.begin(),aux.end());
}

// Min ------------------------------------------
double TeuchosVector::getMinValue( ) const //dummy
{
   const unsigned int size = this->sizeLocal();
   std::vector<double> aux;

  for (unsigned int i=0; i<size; i++ ) {
    aux.push_back((*this)[i]) ;
    }

 return *min_element (aux.begin(),aux.end());

}

// Max index -----------------------------------
int TeuchosVector::getMaxValueIndex( ) const   //dummy
{
  const unsigned int size = this->sizeLocal();
  std::vector<double> vect;

  for (unsigned int i=0; i<size; i++ ) {
    vect.push_back((*this)[i]) ;
  }
  std::vector<double>::iterator iter_max = max_element(vect.begin(), vect.end());
  return distance(vect.begin(), iter_max);
}

// Min index -----------------------------------
int TeuchosVector::getMinValueIndex( ) const //dummy
{
  const unsigned int size = this->sizeLocal();
  std::vector<double> vect;

  for (unsigned int i=0; i<size; i++ ) {
    vect.push_back((*this)[i]) ;
  }
  std::vector<double>::iterator iter_min = min_element(vect.begin(), vect.end());
  return distance(vect.begin(), iter_min);
}

// Max and index -------------------------------
void TeuchosVector::getMaxValueAndIndex( double& max_value, int& max_value_index ) //dummy
{
  const unsigned int size = this->sizeLocal();
  std::vector<double> vect;

  for (unsigned int i=0; i<size; i++ ) {
    vect.push_back((*this)[i]) ;
  }
  std::vector<double>::iterator iter_max = max_element(vect.begin(), vect.end());

  max_value = *iter_max;
  max_value_index =  distance(vect.begin(), iter_max);

  return;
}

// Min and index -------------------------------
void TeuchosVector::getMinValueAndIndex( double& min_value, int& min_value_index ) //dummy
{
  const unsigned int size = this->sizeLocal();
  std::vector<double> vect;

  for (unsigned int i=0; i<size; i++ ) {
    vect.push_back((*this)[i]) ;
  }
  std::vector<double>::iterator iter_min = min_element(vect.begin(), vect.end());

  min_value = *iter_min;
  min_value_index = distance(vect.begin(), iter_min);

  return;
}

// Norm methods ------------------------------------
// -------------------------------------------------
double TeuchosVector::norm2Sq() const
{
 return (m_vec).dot(m_vec );
}

//-------------------------------------------------
double TeuchosVector::norm2() const
{
  return std::sqrt(this->norm2Sq());
}

//-------------------------------------------------
double TeuchosVector::norm1() const
{
  double result = 0.;

  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += fabs((*this)[i]);
  }

  return result;
}

//-------------------------------------------------
double TeuchosVector::normInf() const
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

// Comparison methods -----------------------------
//-------------------------------------------------
bool
TeuchosVector::atLeastOneComponentSmallerThan(const TeuchosVector& rhs) const
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

//-------------------------------------------------
bool
TeuchosVector::atLeastOneComponentBiggerThan(const TeuchosVector& rhs) const
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

//-------------------------------------------------
bool
TeuchosVector::atLeastOneComponentSmallerOrEqualThan(const TeuchosVector& rhs) const
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

//-------------------------------------------------
bool
TeuchosVector::atLeastOneComponentBiggerOrEqualThan(const TeuchosVector& rhs) const
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


// Get/Set methods -----------------------------------
//----------------------------------------------------
void TeuchosVector::cwSet(double value)
{
  (*this)=value;

  return;
}

//----------------------------------------------------
void TeuchosVector::cwSet(unsigned int initialPos, const TeuchosVector& vec)
{
   queso_require_less_msg(initialPos, this->sizeLocal(), "invalid initialPos");

   queso_require_less_equal_msg((initialPos +vec.sizeLocal()), this->sizeLocal(), "invalid vec.sizeLocal()");

  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    (*this)[initialPos+i] = vec[i];
  }

  return;
}


//----------------------------------------------------
/* extracts elements from vector (*this) starting at position initialPos and save in vec */
void TeuchosVector::cwExtract(unsigned int initialPos, TeuchosVector& vec) const
{
   queso_require_less_msg(initialPos, this->sizeLocal(), "invalid initialPos");

   queso_require_less_equal_msg((initialPos +vec.sizeLocal()), this->sizeLocal(), "invalid vec.sizeLocal()");

  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    vec[i] = (*this)[initialPos+i];
  }

  return;
}

//----------------------------------------------------
void TeuchosVector::cwInvert()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = 1./(*this)[i];
  }

  return;
}

//----------------------------------------------------
void TeuchosVector::cwSqrt()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = sqrt((*this)[i]);
  }

  return;
}

//----------------------------------------------------
void
TeuchosVector::cwSetConcatenated(const TeuchosVector& v1, const TeuchosVector& v2)
{
   queso_require_equal_to_msg(this->sizeLocal(), v1.sizeLocal() + v2.sizeLocal(), "incompatible vector sizes");

//  if ( this->sizeLocal() != v1.sizeLocal() + v2.sizeLocal() ) {
//    std::cout << "ERROR in TeuchosVector:: cwSetConcatenated  ---> the vectors' sizes are not compatible.\n";
//    cout << "      in TeuchosVector:: cwSetConcatenated  ---> resizing resulting vector... new size = "
//         << v1.sizeLocal()+v1.sizeLocal() <<endl;
//
//    m_vec.resize(v1.sizeLocal()+v1.sizeLocal());
//  }

  for (unsigned int i = 0; i < v1.sizeLocal(); ++i) {
    (*this)[i] = v1[i];
  }

  for (unsigned int i = 0; i < v2.sizeLocal(); ++i) {
    (*this)[v1.sizeLocal()+i] = v2[i];
  }

  return;
}

// -------------------------------------------------
//updated on 3/18, to use the RngBase+Boost
void TeuchosVector::cwSetGaussian(double mean, double stdDev)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
	(*this)[i] = mean + m_env.rngObject()->gaussianSample(stdDev);
  }
  return;
};

// -------------------------------------------------
//updated on 3/18, to use the RngBase+Boost
void TeuchosVector::cwSetGaussian(const TeuchosVector& meanVec, const TeuchosVector& stdDevVec)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = meanVec[i] + m_env.rngObject()->gaussianSample(stdDevVec[i]);
  }
  return;
};


//----------------------------------------------------
//implemented/checked 2/26/13
//updated on 3/18, to use the RngBase+Boost
 void TeuchosVector::cwSetUniform(const TeuchosVector& aVec, const TeuchosVector& bVec)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = aVec[i] + (bVec[i]-aVec[i])*m_env.rngObject()->uniformSample();
  }
  return;
}


// -------------------------------------------------
//updated on 3/18, to use the RngBase+Boost
void TeuchosVector::cwSetBeta(const TeuchosVector& alpha, const TeuchosVector& beta)
{
  queso_require_equal_to_msg(this->sizeLocal(), alpha.sizeLocal(), "incompatible alpha size");

  queso_require_equal_to_msg(this->sizeLocal(), beta.sizeLocal(), "incompatible beta size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i)
  {
    (*this)[i] = m_env.rngObject()->betaSample(alpha[i],beta[i]);

    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99))
    {
      *m_env.subDisplayFile() << "In TeuchosVector::cwSetBeta()"
                              << ": fullRank "   << m_env.fullRank()
                              << ", i = "        << i
                              << ", alpha[i] = " << alpha[i]
                              << ", beta[i] = "  << beta[i]
                              << ", sample = "   << (*this)[i]
                              << std::endl;
  	}
  }
  return;
};

// -------------------------------------------------
//updated on 3/18, to use the RngBase+Boost
void TeuchosVector::cwSetGamma(const TeuchosVector& aVec, const TeuchosVector& bVec)
{
  queso_require_equal_to_msg(this->sizeLocal(), aVec.sizeLocal(), "incompatible a size");

  queso_require_equal_to_msg(this->sizeLocal(), bVec.sizeLocal(), "incompatible b size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = m_env.rngObject()->gammaSample(aVec[i],bVec[i]);
  }
  return;
}

// -------------------------------------------------
//updated on 3/18, to use the RngBase+Boost
// Using Gamma Distribution to calculate InverseGamma.
// Note the divisions: 1.0/b and the 1.0/generator; they are crucial
void TeuchosVector::cwSetInverseGamma(const TeuchosVector& alpha, const TeuchosVector& beta)
{
  queso_require_equal_to_msg(this->sizeLocal(), alpha.sizeLocal(), "incompatible alpha size");

  queso_require_equal_to_msg(this->sizeLocal(), beta.sizeLocal(), "incompatible beta size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = 1./m_env.rngObject()->gammaSample(alpha[i],1./beta[i]);
  }
  return;
}


// Miscellaneous methods -------------------------
// absolute values -------------------------------
// tested 1/8/13
TeuchosVector
TeuchosVector::abs() const
{
  TeuchosVector abs_of_this_vec( *this );

  unsigned int size = abs_of_this_vec.sizeLocal();

  for( unsigned int i = 0; i < size; ++i )
    {
      abs_of_this_vec[i] = std::fabs( (*this)[i] );
    }

  return abs_of_this_vec;

}

// -------------------------------------------------
void
TeuchosVector::copy_to_std_vector(std::vector<double>& vec)
{
  unsigned int size = this->sizeLocal();
  vec.resize(size);

  for (unsigned int i = 0; i < size; ++i)
  	vec[i] = m_vec[i];

  return;
}

// -------------------------------------------------
void
TeuchosVector::copy_from_std_vector(const std::vector<double> vec)
{
  unsigned int size1 = vec.size(), size2= this->sizeLocal();

  queso_require_equal_to_msg(size1, size2, "vectors have different sizes");

  for (unsigned int i = 0; i < size1; ++i)
      m_vec[i] = vec[i];

  return;
}

// -------------------------------------------------
void TeuchosVector::sort()
{
  std::vector<double> vec;

  (*this).copy_to_std_vector(vec);

  // using default comparison (operator <):
  std::sort (vec.begin(), vec.end());

  (*this).copy_from_std_vector(vec);
};

// -------------------------------------------------
double TeuchosVector::sumOfComponents() const
{
  double result = 0.;
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += (*this)[i];
  }

  return result;
}

// -------------------------------------------------
// added 2/28/13
void
TeuchosVector::mpiBcast(int srcRank, const MpiComm& bcastComm)
{
  // Filter out those nodes that should not participate
  if (bcastComm.MyPID() < 0) return;

  // Check 'srcRank'
  queso_require_msg(!((srcRank < 0) || (srcRank >= bcastComm.NumProc())), "invalud srcRank");

  // Check number of participant nodes
  double localNumNodes = 1.;
  double totalNumNodes = 0.;
  bcastComm.Allreduce<double>(&localNumNodes, &totalNumNodes, (int) 1, RawValue_MPI_SUM,
                      "TeuchosVector::mpiBcast()",
                      "failed MPI.Allreduce() for numNodes");
  queso_require_equal_to_msg(((int) totalNumNodes), bcastComm.NumProc(), "inconsistent numNodes");

  // Check that all participant nodes have the same vector size
  double localVectorSize  = this->sizeLocal();
  double sumOfVectorSizes = 0.;
  bcastComm.Allreduce<double>(&localVectorSize, &sumOfVectorSizes, (int) 1, RawValue_MPI_SUM,
                      "TeuchosVector::mpiBcast()",
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
                  "TeuchosVector::mpiBcast()",
                  "failed MPI.Bcast()");

  if (bcastComm.MyPID() != srcRank) {
    for (unsigned int i = 0; i < dataBuffer.size(); ++i) {
      (*this)[i] = dataBuffer[i];
    }
  }

  return;
}

// -------------------------------------------------
// added 2/28/13
void
TeuchosVector::mpiAllReduce(RawType_MPI_Op mpiOperation, const MpiComm& opComm, TeuchosVector& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  unsigned int size = this->sizeLocal();
  queso_require_equal_to_msg(size, resultVec.sizeLocal(), "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double srcValue = (*this)[i];
    double resultValue = 0.;
    opComm.Allreduce<double>(&srcValue, &resultValue, (int) 1, mpiOperation,
                     "TeuchosVector::mpiAllReduce()",
                     "failed MPI.Allreduce()");
    resultVec[i] = resultValue;
  }

  return;
}

// -------------------------------------------------
// added 2/28/13
void
TeuchosVector::mpiAllQuantile(double probability, const MpiComm& opComm, TeuchosVector& resultVec) const
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
                  "TeuchosVector::mpiAllQuantile()",
                  "failed MPI.Gather()");

    std::sort(vecOfDoubles.begin(), vecOfDoubles.end());

    double result = vecOfDoubles[(unsigned int)( probability*((double)(vecOfDoubles.size()-1)) )];

    opComm.Bcast((void *) &result, (int) 1, RawValue_MPI_DOUBLE, 0,
                 "TeuchosVector::mpiAllQuantile()",
                 "failed MPI.Bcast()");

    resultVec[i] = result;
  }

  return;
}

// -------------------------------------------------
// added/tested 2/28/13
void
TeuchosVector::matlabLinearInterpExtrap(
  const TeuchosVector& x1Vec,
  const TeuchosVector& y1Vec,
  const TeuchosVector& x2Vec)
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

// -------------------------------------------------
// added 2/28/13
void
TeuchosVector::matlabDiff(
  unsigned int      firstPositionToStoreDiff,
  double            valueForRemainderPosition,
  TeuchosVector& outputVec) const
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

// I/O methods -------------------------------------
// -------------------------------------------------
void
TeuchosVector::print(std::ostream& os) const
{
  unsigned int size = this->sizeLocal();

  std::ostream::fmtflags curr_fmt = os.flags();

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

  os.flags(curr_fmt);
  return;
}

// -------------------------------------------------
void
TeuchosVector::subReadContents(
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
    queso_require_equal_to_msg(tmpString, "=", "string should be the '=' sign");

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
      *m_env.subDisplayFile() << "In TeuchosVector::subReadContents()"
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
      *m_env.subDisplayFile() << "In TeuchosVector::subReadContents()"
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
    queso_require_equal_to_msg(tmpString, "=", "in core 0, string should be the '=' sign");

    // Take into account the ' [' portion
    std::streampos tmpPos = filePtrSet.ifsVar->tellg();
    filePtrSet.ifsVar->seekg(tmpPos+(std::streampos)2);

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In TeuchosVector::subReadContents()"
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

// -------------------------------------------------
void
TeuchosVector::subWriteContents(
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


// -------------------------------------------------
// Private metfods ---------------------------------
// -------------------------------------------------

void
TeuchosVector::copy(const TeuchosVector& rhs)
{
  this->Vector::base_copy(rhs);

  unsigned int size1 = m_vec.length();
  unsigned int size2 = rhs.sizeLocal();

  if (size1==size2){
    for (unsigned int i=0;i<size1;i++){
    m_vec[i]=rhs[i];
    }
  }
  return;
}


// -------------------------------------------------
// Operators outside class definition --------------
// -------------------------------------------------

// -------------------------------------------------
TeuchosVector operator/(double a, const TeuchosVector& x)
{
  TeuchosVector answer(x); //copy x to answer
  answer.cwInvert();
  answer *= a;
  return answer;
}

// -------------------------------------------------
TeuchosVector operator/(const TeuchosVector& x, const TeuchosVector& y)
{
  TeuchosVector answer(x);
  answer /= y;
  return answer;
}

// -------------------------------------------------
TeuchosVector operator*(double a, const TeuchosVector& x)
{
  TeuchosVector answer(x);
  answer *= a;
  return answer;
}

// -------------------------------------------------
TeuchosVector operator*(const TeuchosVector& x, const TeuchosVector& y)
{
  TeuchosVector answer(x);
  answer *= y;
  return answer;
}

// -------------------------------------------------
double scalarProduct(const TeuchosVector& x, const TeuchosVector& y)
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

// -------------------------------------------------
TeuchosVector operator+(const TeuchosVector& x, const TeuchosVector& y)
{
  TeuchosVector answer(x);
  answer += y;
  return answer;
}

// -------------------------------------------------
TeuchosVector operator-(const TeuchosVector& x, const TeuchosVector& y)
{
  TeuchosVector answer(x);
  answer -= y;
  return answer;
}

// -------------------------------------------------
bool
operator== (const TeuchosVector& lhs, const TeuchosVector& rhs)
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

// -------------------------------------------------
std::ostream&
operator<<(std::ostream& os, const TeuchosVector& obj)
{
  obj.print(os);

  return os;
}

// -------------------------------------------------

}  // End namespace QUESO

#endif // ifdef QUESO_HAS_TRILINOS
