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

#include <uqTeuchosVector.h>

namespace QUESO {

using std:: cout;
using std:: endl;

// ---------------------------------------------------
// default constructor ------------------------------- 
TeuchosVectorClass::TeuchosVectorClass()  :
  VectorClass()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(), default",
                      "should not be used by user");  
};


// constructor with dimension ----------------------- 
TeuchosVectorClass::TeuchosVectorClass(const BaseEnvironmentClass& env, const MapClass& map)
  :
  VectorClass(env,map)
{
  m_vec.size(map.NumGlobalElements());

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(1)",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() != map.NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(1)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() != map.NumGlobalElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(1)",
                      "incompatible global vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() != m_map.NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(1)",
                      "incompatible own vec size");

  //std::cout << "In TeuchosVectorClass::constructor(env,map)"
  //          << "\n  m_vec.length()             = " << m_vec.length()
  //          << "\n  map.NumGlobalElements() = " << map.NumGlobalElements()
  //          << "\n  map.NumMyElements()     = " << map.NumMyElements()
  //          << std::endl;
}

// ---------------------------------------------------
TeuchosVectorClass::TeuchosVectorClass(const BaseEnvironmentClass& env, const MapClass& map, double value)
  :
  VectorClass(env,map)
{
  m_vec.size(map.NumGlobalElements());
  m_vec = value;
    
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(2)",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  map.NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(2)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  map.NumGlobalElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(2)",
                      "incompatible global vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  m_map.NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(2)",
                      "incompatible own vec size");
}

// ---------------------------------------------------
TeuchosVectorClass::TeuchosVectorClass(const BaseEnvironmentClass& env, double d1, double d2, const MapClass& map)
  :
  VectorClass(env,map)
  {
  m_vec.size(map.NumGlobalElements());

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(3), linspace",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  map.NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(3)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  map.NumGlobalElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(3)",
                      "incompatible global vec size");

  for (int i = 0; i < m_vec.length(); ++i) {
    double alpha = (double) i / ((double) m_vec.length() - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  m_map.NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(3)",
                      "incompatible own vec size");
}

// ---------------------------------------------------
TeuchosVectorClass::TeuchosVectorClass(const TeuchosVectorClass& v, double d1, double d2)
  :
  VectorClass(v.env(),v.map())
{
  m_vec.size(v.sizeLocal());
  
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(4), linspace",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  v.map().NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(4)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  v.map().NumGlobalElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(4)",
                      "incompatible global vec size");

  for ( int i = 0; i < m_vec.length(); ++i) {
    double alpha = (double) i / ((double) m_vec.length() - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  m_map.NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(4)",
                      "incompatible own vec size");
}

// ---------------------------------------------------
TeuchosVectorClass::TeuchosVectorClass(const TeuchosVectorClass& v)  // mox
  :
  VectorClass(v.env(),v.map())
 {
   m_vec.size(v.sizeLocal());
   
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(5), copy",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  v.map().NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(5)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  v.map().NumGlobalElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(5)",
                      "incompatible global vec size");
  this->copy(v);

  UQ_FATAL_TEST_MACRO(m_vec.length() != m_map.NumMyElements(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::constructor(5)",
                      "incompatible own vec size");
}


// destructor -------------------------------------- 
TeuchosVectorClass::~TeuchosVectorClass()
{
};

// Set methods ------------------------------------
//-------------------------------------------------
// assigns values a to all entrances in the vector
TeuchosVectorClass& TeuchosVectorClass::operator=(double a)
{
  m_vec.putScalar(a);
  return *this;
}

//-------------------------------------------------
TeuchosVectorClass& TeuchosVectorClass::operator=(const TeuchosVectorClass& rhs)
{
  unsigned int size1 = m_vec.length();
  unsigned int size2 = rhs.sizeLocal();
  
  UQ_FATAL_TEST_MACRO(size1 != size2, // mox
                      m_env.worldRank(),
                      "TeuchosVectorClass::operator=()",
                      "the vectors do NOT have the same size.\n");
               
  if (size1==size2){
    for (unsigned int i=0;i<size1;i++){
    m_vec[i]=rhs[i];  
    }
  }

  return *this;
}

//-------------------------------------------------
TeuchosVectorClass& TeuchosVectorClass::operator*=(double a)
{
  m_vec.scale(a); 	//Scale this vector by a; *this = a*this. 
  return *this;
}

 //-------------------------------------------------
TeuchosVectorClass& TeuchosVectorClass::operator/=(double a)
{
  (*this) *= (1.0/a);
  return *this;
}

//-------------------------------------------------
TeuchosVectorClass& TeuchosVectorClass::operator*=(const TeuchosVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      m_env.worldRank(),
                      "TeuchosVectorClass::operator*=()",
                      "the vectors do NOT have the same size.\n");
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] *= rhs[i];
    }
 }
 return *this;
}

//-------------------------------------------------
TeuchosVectorClass& TeuchosVectorClass::operator/=(const TeuchosVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      m_env.worldRank(),
                      "TeuchosVectorClass::operator/=()",
                      "the vectors do NOT have the same size.\n");
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] /= rhs[i];
    }
  }
  return *this;
}

//-------------------------------------------------
TeuchosVectorClass& TeuchosVectorClass::operator+=(const TeuchosVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();

  UQ_FATAL_TEST_MACRO((size1 != size2),
                      m_env.worldRank(),
                      "TeuchosVectorClass::operator+=()",
                      "the vectors do NOT have the same size.\n");
  
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] += rhs[i];
    }
  }
  return *this;
}

//-------------------------------------------------
TeuchosVectorClass& TeuchosVectorClass::operator-=(const TeuchosVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();

  UQ_FATAL_TEST_MACRO((size1 != size2),
                      m_env.worldRank(),
                      "TeuchosVectorClass::operator-=()",
                      "the vectors do NOT have the same size.\n");
  
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] -= rhs[i];
    }
  }
  
  return *this;
}


// Accessor methods --------------------------------
//-------------------------------------------------
double& TeuchosVectorClass::operator[](unsigned int i)
{
  return m_vec[i];
}

//-------------------------------------------------
const double& TeuchosVectorClass::operator[](unsigned int i) const
{
  return m_vec[i];
}

// Attribute methods ------------------------------
//-------------------------------------------------
unsigned int TeuchosVectorClass::sizeLocal() const
{
  UQ_FATAL_TEST_MACRO(m_vec.length() != (int) m_map.NumMyElements(),
                       m_env.worldRank(),
                       "TeuchosVectorClass::sizeLocal()",
                       "incompatible vec size");
  
  return m_vec.length();
}

//-------------------------------------------------
unsigned int TeuchosVectorClass::sizeGlobal() const
{
   UQ_FATAL_TEST_MACRO(m_vec.length() != (int) m_map.NumGlobalElements(),
                       m_env.worldRank(),
                       "TeuchosVectorClass::sizeGlobal()",
                       "incompatible vec size");
  return m_vec.length();
}

//-------------------------------------------------
// TODO: needs to be checked. It may not be used at all. Kemelli 4/30/13.
double*
TeuchosVectorClass::values()
{
  return  m_vec.values();
};

// Getting Max and Min values -------------------
// Max ------------------------------------------
double TeuchosVectorClass::getMaxValue( ) const  //dummy
{
  const unsigned int size = this->sizeLocal();
  std::vector<double> aux;
    
  for (unsigned int i=0; i<size; i++ ) {
    aux.push_back((*this)[i]) ;
  }
 
  return *max_element (aux.begin(),aux.end());
}

// Min ------------------------------------------
double TeuchosVectorClass::getMinValue( ) const //dummy
{  
   const unsigned int size = this->sizeLocal();
   std::vector<double> aux;
  
  for (unsigned int i=0; i<size; i++ ) {
    aux.push_back((*this)[i]) ;
    }
  
 return *min_element (aux.begin(),aux.end());
  
}

// Max index -----------------------------------
int TeuchosVectorClass::getMaxValueIndex( ) const   //dummy
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
int TeuchosVectorClass::getMinValueIndex( ) const //dummy
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
void TeuchosVectorClass::getMaxValueAndIndex( double& max_value, int& max_value_index ) //dummy
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
void TeuchosVectorClass::getMinValueAndIndex( double& min_value, int& min_value_index ) //dummy
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
double TeuchosVectorClass::norm2Sq() const
{  
 return (m_vec).dot(m_vec );  
}

//-------------------------------------------------
double TeuchosVectorClass::norm2() const
{
  return std::sqrt(this->norm2Sq());
}

//-------------------------------------------------
double TeuchosVectorClass::norm1() const
{
  double result = 0.;

  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += fabs((*this)[i]);
  }

  return result;
}

//-------------------------------------------------
double TeuchosVectorClass::normInf() const
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
TeuchosVectorClass::atLeastOneComponentSmallerThan(const TeuchosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "TeuchosVectorClass::atLeastOneComponentSmallerThan()",
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

//-------------------------------------------------
bool
TeuchosVectorClass::atLeastOneComponentBiggerThan(const TeuchosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "TeuchosVectorClass::atLeastOneComponentBiggerThan()",
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

//-------------------------------------------------
bool
TeuchosVectorClass::atLeastOneComponentSmallerOrEqualThan(const TeuchosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "TeuchosVectorClass::atLeastOneComponentSmallerOrEqualThan()",
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

//-------------------------------------------------
bool
TeuchosVectorClass::atLeastOneComponentBiggerOrEqualThan(const TeuchosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "TeuchosVectorClass::atLeastOneComponentBiggerOrEqualThan()",
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


// Get/Set methods -----------------------------------
//----------------------------------------------------
void TeuchosVectorClass::cwSet(double value)
{
  (*this)=value;
  
  return;
}

//----------------------------------------------------
void TeuchosVectorClass::cwSet(unsigned int initialPos, const TeuchosVectorClass& vec)
{
   UQ_FATAL_TEST_MACRO(initialPos >= this->sizeLocal(),
                       m_env.worldRank(),
                       "TeuchosVectorClass::cwSet()",
                       "invalid initialPos");
 
   UQ_FATAL_TEST_MACRO((initialPos +vec.sizeLocal()) > this->sizeLocal(),
                       m_env.worldRank(),
                       "TeuchosVectorClass::cwSet()",
                       "invalid vec.sizeLocal()");
  
  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    (*this)[initialPos+i] = vec[i];
  }

  return;
}


//----------------------------------------------------
/* extracts elements from vector (*this) starting at position initialPos and save in vec */
void TeuchosVectorClass::cwExtract(unsigned int initialPos, TeuchosVectorClass& vec) const
{
   UQ_FATAL_TEST_MACRO(initialPos >= this->sizeLocal(),
                       m_env.worldRank(),
                       "TeuchosVectorClass::cwExtract()",
                       "invalid initialPos");
 
   UQ_FATAL_TEST_MACRO((initialPos +vec.sizeLocal()) > this->sizeLocal(),
                       m_env.worldRank(),
                       "TeuchosVectorClass::cwExtract()",
                       "invalid vec.sizeLocal()");
 
  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    vec[i] = (*this)[initialPos+i];
  }

  return;
}

//----------------------------------------------------
void TeuchosVectorClass::cwInvert()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = 1./(*this)[i];
  }

  return;
}

//----------------------------------------------------
void TeuchosVectorClass::cwSqrt()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = sqrt((*this)[i]);
  }

  return;
}

//----------------------------------------------------
void
TeuchosVectorClass::cwSetConcatenated(const TeuchosVectorClass& v1, const TeuchosVectorClass& v2)
{
   UQ_FATAL_TEST_MACRO(this->sizeLocal() != v1.sizeLocal() + v2.sizeLocal(),
                       m_env.worldRank(),
                       "TeuchosVectorClass::cwSetConcatenated(1)",
                       "incompatible vector sizes");

//  if ( this->sizeLocal() != v1.sizeLocal() + v2.sizeLocal() ) {
//    std::cout << "ERROR in TeuchosVectorClass:: cwSetConcatenated  ---> the vectors' sizes are not compatible.\n";
//    cout << "      in TeuchosVectorClass:: cwSetConcatenated  ---> resizing resulting vector... new size = " 
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
//updated on 3/18, to use the RngBaseClass+Boost
void TeuchosVectorClass::cwSetGaussian(double mean, double stdDev)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
	(*this)[i] = mean + m_env.rngObject()->gaussianSample(stdDev);
  }
  return;
};

// -------------------------------------------------
//updated on 3/18, to use the RngBaseClass+Boost
void TeuchosVectorClass::cwSetGaussian(const TeuchosVectorClass& meanVec, const TeuchosVectorClass& stdDevVec)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = meanVec[i] + m_env.rngObject()->gaussianSample(stdDevVec[i]);
  }  
  return;
};


//----------------------------------------------------
//implemented/checked 2/26/13
//updated on 3/18, to use the RngBaseClass+Boost
 void TeuchosVectorClass::cwSetUniform(const TeuchosVectorClass& aVec, const TeuchosVectorClass& bVec)
{
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = aVec[i] + (bVec[i]-aVec[i])*m_env.rngObject()->uniformSample();
  }
  return;
}


// -------------------------------------------------
//updated on 3/18, to use the RngBaseClass+Boost
void TeuchosVectorClass::cwSetBeta(const TeuchosVectorClass& alpha, const TeuchosVectorClass& beta)
{
  UQ_FATAL_TEST_MACRO(this->sizeLocal() != alpha.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::cwSetBeta()",
                      "incompatible alpha size");

  UQ_FATAL_TEST_MACRO(this->sizeLocal() != beta.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::cwSetBeta()",
                      "incompatible beta size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) 
  {
    (*this)[i] = m_env.rngObject()->betaSample(alpha[i],beta[i]);
    
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) 
    {
      *m_env.subDisplayFile() << "In TeuchosVectorClass::cwSetBeta()"
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
//updated on 3/18, to use the RngBaseClass+Boost
void TeuchosVectorClass::cwSetGamma(const TeuchosVectorClass& aVec, const TeuchosVectorClass& bVec)
{
  UQ_FATAL_TEST_MACRO(this->sizeLocal() != aVec.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::cwSetGamma()",
                      "incompatible a size");

  UQ_FATAL_TEST_MACRO(this->sizeLocal() != bVec.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::cwSetGamma()",
                      "incompatible b size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = m_env.rngObject()->gammaSample(aVec[i],bVec[i]);
  }
  return;
}

// -------------------------------------------------
//updated on 3/18, to use the RngBaseClass+Boost
// Using Gamma Distribution to calculate InverseGamma.
// Note the divisions: 1.0/b and the 1.0/generator; they are crucial
void TeuchosVectorClass::cwSetInverseGamma(const TeuchosVectorClass& alpha, const TeuchosVectorClass& beta)
{
  UQ_FATAL_TEST_MACRO(this->sizeLocal() != alpha.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::cwSetInverseGamma()",
                      "incompatible alpha size");

  UQ_FATAL_TEST_MACRO(this->sizeLocal() != beta.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::cwSetInverseGamma()",
                      "incompatible beta size");

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) {
    (*this)[i] = 1./m_env.rngObject()->gammaSample(alpha[i],1./beta[i]);
  }
  return;
}


// Miscellaneous methods -------------------------
// absolute values -------------------------------
// tested 1/8/13
TeuchosVectorClass
TeuchosVectorClass::abs() const
{
  TeuchosVectorClass abs_of_this_vec( *this );

  unsigned int size = abs_of_this_vec.sizeLocal();

  for( unsigned int i = 0; i < size; ++i )
    {
      abs_of_this_vec[i] = std::fabs( (*this)[i] );
    }

  return abs_of_this_vec;

}

// -------------------------------------------------
void
TeuchosVectorClass::copy_to_std_vector(std::vector<double>& vec)
{
  unsigned int size = this->sizeLocal();
  vec.resize(size);
  
  for (unsigned int i = 0; i < size; ++i) 
  	vec[i] = m_vec[i];
    
  return;
}

// -------------------------------------------------
void
TeuchosVectorClass::copy_from_std_vector(const std::vector<double> vec)
{
  unsigned int size1 = vec.size(), size2= this->sizeLocal();
 
  UQ_FATAL_TEST_MACRO((size1 != size2),
                       m_env.worldRank(),
                       "In TeuchosVectorClass::copy_from_std_vector()",
					   "vectors have different sizes");
 
  for (unsigned int i = 0; i < size1; ++i) 
      m_vec[i] = vec[i];
  
  return;
}

// -------------------------------------------------
void TeuchosVectorClass::sort()
{
  std::vector<double> vec;
  
  (*this).copy_to_std_vector(vec);  
  
  // using default comparison (operator <):
  std::sort (vec.begin(), vec.end()); 
  
  (*this).copy_from_std_vector(vec); 
};

// -------------------------------------------------
double TeuchosVectorClass::sumOfComponents() const
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
TeuchosVectorClass::mpiBcast(int srcRank, const MpiCommClass& bcastComm)
{
  // Filter out those nodes that should not participate
  if (bcastComm.MyPID() < 0) return;

  // Check 'srcRank'
  UQ_FATAL_TEST_MACRO((srcRank < 0) || (srcRank >= bcastComm.NumProc()),
                      m_env.worldRank(),
                      "TeuchosVectorClass::mpiBcast()",
                      "invalud srcRank");

  // Check number of participant nodes
  double localNumNodes = 1.;
  double totalNumNodes = 0.;
  bcastComm.Allreduce((void *) &localNumNodes, (void *) &totalNumNodes, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_SUM,
                      "TeuchosVectorClass::mpiBcast()",
                      "failed MPI.Allreduce() for numNodes");
  UQ_FATAL_TEST_MACRO(((int) totalNumNodes) != bcastComm.NumProc(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::mpiBcast()",
                      "inconsistent numNodes");

  // Check that all participant nodes have the same vector size
  double localVectorSize  = this->sizeLocal();
  double sumOfVectorSizes = 0.; 
  bcastComm.Allreduce((void *) &localVectorSize, (void *) &sumOfVectorSizes, (int) 1, RawValue_MPI_DOUBLE, RawValue_MPI_SUM,
                      "TeuchosVectorClass::mpiBcast()",
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
                      "TeuchosVectorClass::mpiBcast()",
                      "inconsistent vectorSize");

  // Ok, bcast data
  std::vector<double> dataBuffer((unsigned int) localVectorSize, 0.);
  if (bcastComm.MyPID() == srcRank) {
    for (unsigned int i = 0; i < dataBuffer.size(); ++i) {
      dataBuffer[i] = (*this)[i];
    }
  }

  bcastComm.Bcast((void *) &dataBuffer[0], (int) localVectorSize, RawValue_MPI_DOUBLE, srcRank,
                  "TeuchosVectorClass::mpiBcast()",
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
TeuchosVectorClass::mpiAllReduce(RawType_MPI_Op mpiOperation, const MpiCommClass& opComm, TeuchosVectorClass& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  unsigned int size = this->sizeLocal();
  UQ_FATAL_TEST_MACRO(size != resultVec.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::mpiAllReduce()",
                      "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double srcValue = (*this)[i];
    double resultValue = 0.;
    opComm.Allreduce((void *) &srcValue, (void *) &resultValue, (int) 1, RawValue_MPI_DOUBLE, mpiOperation,
                     "TeuchosVectorClass::mpiAllReduce()",
                     "failed MPI.Allreduce()");
    resultVec[i] = resultValue;
  }

  return;
}

// -------------------------------------------------
// added 2/28/13
void
TeuchosVectorClass::mpiAllQuantile(double probability, const MpiCommClass& opComm, TeuchosVectorClass& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  UQ_FATAL_TEST_MACRO((probability < 0.) || (1. < probability),
                      m_env.worldRank(),
                      "TeuchosVectorClass::mpiAllQuantile()",
                      "invalid input");

  unsigned int size = this->sizeLocal();
  UQ_FATAL_TEST_MACRO(size != resultVec.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::mpiAllQuantile()",
                      "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double auxDouble = (int) (*this)[i];
    std::vector<double> vecOfDoubles(opComm.NumProc(),0.);
    opComm.Gather((void *) &auxDouble, 1, RawValue_MPI_DOUBLE, (void *) &vecOfDoubles[0], (int) 1, RawValue_MPI_DOUBLE, 0,
                  "TeuchosVectorClass::mpiAllQuantile()",
                  "failed MPI.Gather()");

    std::sort(vecOfDoubles.begin(), vecOfDoubles.end());

    double result = vecOfDoubles[(unsigned int)( probability*((double)(vecOfDoubles.size()-1)) )];

    opComm.Bcast((void *) &result, (int) 1, RawValue_MPI_DOUBLE, 0,
                 "TeuchosVectorClass::mpiAllQuantile()",
                 "failed MPI.Bcast()");

    resultVec[i] = result;
  }

  return;
}

// -------------------------------------------------
// added/tested 2/28/13
void
TeuchosVectorClass::matlabLinearInterpExtrap(
  const TeuchosVectorClass& x1Vec,
  const TeuchosVectorClass& y1Vec,
  const TeuchosVectorClass& x2Vec)
{
  UQ_FATAL_TEST_MACRO(x1Vec.sizeLocal() <= 1,
                      m_env.worldRank(),
                      "TeuchosVectorClass::matlabLinearInterpExtrap()",
                      "invalid 'x1' size");

  UQ_FATAL_TEST_MACRO(x1Vec.sizeLocal() != y1Vec.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::matlabLinearInterpExtrap()",
                      "invalid 'x1' and 'y1' sizes");

  UQ_FATAL_TEST_MACRO(x2Vec.sizeLocal() != this->sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::matlabLinearInterpExtrap()",
                      "invalid 'x2' and 'this' sizes");

  for (unsigned int i = 1; i < x1Vec.sizeLocal(); ++i) { // Yes, '1'
    UQ_FATAL_TEST_MACRO(x1Vec[i] <= x1Vec[i-1],
                        m_env.worldRank(),
                        "TeuchosVectorClass::matlabLinearInterpExtrap()",
                        "invalid 'x1' values");
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
TeuchosVectorClass::matlabDiff(
  unsigned int      firstPositionToStoreDiff,
  double            valueForRemainderPosition,
  TeuchosVectorClass& outputVec) const
{
  unsigned int size = this->sizeLocal();

  UQ_FATAL_TEST_MACRO(firstPositionToStoreDiff > 1,
                      m_env.worldRank(),
                      "TeuchosVectorClass::matlabDiff()",
                      "invalid firstPositionToStoreDiff");

  UQ_FATAL_TEST_MACRO(size != outputVec.sizeLocal(),
                      m_env.worldRank(),
                      "TeuchosVectorClass::matlabDiff()",
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

// I/O methods -------------------------------------
// -------------------------------------------------
void
TeuchosVectorClass::print(std::ostream& os) const
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
TeuchosVectorClass::subReadContents(
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds)
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "TeuchosVectorClass::subReadContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "TeuchosVectorClass::subReadContents()",
                      "implemented just for sequential vectors for now");

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
    UQ_FATAL_TEST_MACRO(tmpString != "=",
                        m_env.worldRank(),
                        "TeuchosVectorClass::subReadContents()",
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
                          "TeuchosVectorClass::subReadContents()",
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
                          "TeuchosVectorClass::subReadContents()",
                          "symbol ')' not found in first line of file");
      nParamsString[posInParamsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ')');
    nParamsString[posInParamsString] = '\0';

    // Convert 'n_positions' and 'n_params' strings to numbers
    unsigned int sizeOfVecInFile = (unsigned int) strtod(nPositionsString,NULL);
    unsigned int numParamsInFile = (unsigned int) strtod(nParamsString,   NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In TeuchosVectorClass::subReadContents()"
                              << ": fullRank "            << m_env.fullRank()
                              << ", sizeOfVecInFile = "   << sizeOfVecInFile
                              << ", numParamsInFile = "   << numParamsInFile
                              << ", this->sizeLocal() = " << this->sizeLocal()
                              << std::endl;
    }

    // Check if [size of vec in file] >= [requested sub vec size]
    UQ_FATAL_TEST_MACRO(sizeOfVecInFile < subReadSize,
                        m_env.worldRank(),
                        "TeuchosVectorClass::subReadContents()",
                        "size of vec in file is not big enough");

    // Check if [num params in file] == [num params in current vec]
    UQ_FATAL_TEST_MACRO(numParamsInFile != numParams,
                        m_env.worldRank(),
                        "TeuchosVectorClass::subReadContents()",
                        "number of parameters of vec in file is different than number of parameters in this vec object");

    // Code common to any core in a communicator
    unsigned int maxCharsPerLine = 64*numParams; // Up to about 60 characters to represent each parameter value

    unsigned int lineId = 0;
    while (lineId < idOfMyFirstLine) {
      filePtrSet.ifsVar->ignore(maxCharsPerLine,'\n');
      lineId++;
    };

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In TeuchosVectorClass::subReadContents()"
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
                        "TeuchosVectorClass::subReadContents()",
                        "in core 0, string should be the '=' sign");

    // Take into account the ' [' portion
    std::streampos tmpPos = filePtrSet.ifsVar->tellg();
    filePtrSet.ifsVar->seekg(tmpPos+(std::streampos)2);

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In TeuchosVectorClass::subReadContents()"
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
TeuchosVectorClass::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "TeuchosVectorClass::subWriteContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "TeuchosVectorClass::subWriteContents()",
                      "implemented just for sequential vectors for now");

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
TeuchosVectorClass::copy(const TeuchosVectorClass& rhs)
{
  this->VectorClass::copy(rhs); 
  
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
TeuchosVectorClass operator/(double a, const TeuchosVectorClass& x)
{
  TeuchosVectorClass answer(x); //copy x to answer
  answer.cwInvert();
  answer *= a;
  return answer;
}

// -------------------------------------------------
TeuchosVectorClass operator/(const TeuchosVectorClass& x, const TeuchosVectorClass& y)
{
  TeuchosVectorClass answer(x);
  answer /= y;
  return answer;
}

// -------------------------------------------------
TeuchosVectorClass operator*(double a, const TeuchosVectorClass& x)
{
  TeuchosVectorClass answer(x);
  answer *= a;
  return answer;
}

// -------------------------------------------------
TeuchosVectorClass operator*(const TeuchosVectorClass& x, const TeuchosVectorClass& y)
{
  TeuchosVectorClass answer(x);
  answer *= y;
  return answer;
}

// -------------------------------------------------
double scalarProduct(const TeuchosVectorClass& x, const TeuchosVectorClass& y)
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

// -------------------------------------------------
TeuchosVectorClass operator+(const TeuchosVectorClass& x, const TeuchosVectorClass& y)
{
  TeuchosVectorClass answer(x);
  answer += y;
  return answer;
}

// -------------------------------------------------
TeuchosVectorClass operator-(const TeuchosVectorClass& x, const TeuchosVectorClass& y)
{
  TeuchosVectorClass answer(x);
  answer -= y;
  return answer;
}

// -------------------------------------------------
bool
operator== (const TeuchosVectorClass& lhs, const TeuchosVectorClass& rhs)
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

// -------------------------------------------------
std::ostream&
operator<<(std::ostream& os, const TeuchosVectorClass& obj)
{
  obj.print(os);

  return os;
}

// -------------------------------------------------

}  // End namespace QUESO

#endif // ifdef QUESO_HAS_TRILINOS
