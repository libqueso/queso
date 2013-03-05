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

#include <uqDefines.h>

#ifdef QUESO_HAS_TRILINOS

#include <uqTeuchosVector.h>

using std:: cout;
using std:: endl;

// standard constructor ----------------------------- 
uqTeuchosVectorClass::uqTeuchosVectorClass()  :
  uqVectorClass()
{
  //UQ_FATAL_TEST_MACRO(true,
  //                    m_env.worldRank(),
  //                    "uqTeuchosVectorClass::constructor(), default",
  //                    "should not be used by user");
  cout << "Constructor 1" << endl;
};


// constructor with dimension ----------------------- 

uqTeuchosVectorClass::uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map)
  :
  uqVectorClass(env,map)
{
  m_vec.size(map.NumGlobalElements());
  //std::cout << "Entering uqTeuchosVectorClass::constructor(1)" << std::endl;

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(1)",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() != map.NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(1)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() != map.NumGlobalElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(1)",
                      "incompatible global vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() != m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(1)",
                      "incompatible own vec size");

  //std::cout << "In uqTeuchosVectorClass::constructor(env,map)"
  //          << "\n  m_vec.length()             = " << m_vec.length()
  //          << "\n  map.NumGlobalElements() = " << map.NumGlobalElements()
  //          << "\n  map.NumMyElements()     = " << map.NumMyElements()
  //          << std::endl;

  //std::cout << "Leaving uqTeuchosVectorClass::constructor(1)" << std::endl;
}

uqTeuchosVectorClass::uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double value)
  :
  uqVectorClass(env,map)
 
{
  m_vec.size(map.NumGlobalElements());
  m_vec = value;
    
  //std::cout << "Entering uqTeuchosVectorClass::constructor(2)" << std::endl;

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(2)",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  map.NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(2)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  map.NumGlobalElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(2)",
                      "incompatible global vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(2)",
                      "incompatible own vec size");

  //std::cout << "Leaving uqTeuchosVectorClass::constructor(2)" << std::endl;
}


uqTeuchosVectorClass::uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const uqMapClass& map)
  :
  uqVectorClass(env,map)
  {
  m_vec.size(map.NumGlobalElements());

  //std::cout << "Entering uqTeuchosVectorClass::constructor(3)" << std::endl;

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(3), linspace",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  map.NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(3)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  map.NumGlobalElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(3)",
                      "incompatible global vec size");

  for (int i = 0; i < m_vec.length(); ++i) {
    double alpha = (double) i / ((double) m_vec.length() - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(3)",
                      "incompatible own vec size");

  //std::cout << "Leaving uqTeuchosVectorClass::constructor(3)" << std::endl;
}


uqTeuchosVectorClass::uqTeuchosVectorClass(const uqTeuchosVectorClass& v, double d1, double d2)
  :
  uqVectorClass(v.env(),v.map())
{
  m_vec.size(v.sizeLocal());
  
  //std::cout << "Entering uqTeuchosVectorClass::constructor(4)" << std::endl;

  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(4), linspace",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  v.map().NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(4)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  v.map().NumGlobalElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(4)",
                      "incompatible global vec size");

  for ( int i = 0; i < m_vec.length(); ++i) {
    double alpha = (double) i / ((double) m_vec.length() - 1.);
    (*this)[i] = (1.-alpha)*d1 + alpha*d2;
  }

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(4)",
                      "incompatible own vec size");

  //std::cout << "Leaving uqTeuchosVectorClass::constructor(4)" << std::endl;
}



uqTeuchosVectorClass::uqTeuchosVectorClass(const uqTeuchosVectorClass& v)  // mox
  :
  uqVectorClass(v.env(),v.map())
 {
   m_vec.size(v.sizeLocal());
  //std::cout << "Entering uqTeuchosVectorClass::constructor(5)" << std::endl;

  // prudenci 2010-06-17 mox
  UQ_FATAL_TEST_MACRO((m_vec == NULL),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(5), copy",
                      "null vector generated");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  v.map().NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(5)",
                      "incompatible local vec size");

  UQ_FATAL_TEST_MACRO(m_vec.length() !=  v.map().NumGlobalElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(5)",
                      "incompatible global vec size");

  this->copy(v);

  UQ_FATAL_TEST_MACRO(m_vec.length() != m_map.NumMyElements(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::constructor(5)",
                      "incompatible own vec size");

  //std::cout << "Leaving uqTeuchosVectorClass::constructor(5)" << std::endl;
}


// destructor -------------------------------------- 
uqTeuchosVectorClass::~uqTeuchosVectorClass()
{
};


// -------------------------------------------------
//TODO : find a smart way to define seed
void uqTeuchosVectorClass::cwSetGaussian(double mean, double stdDev)
{
  int seed =1;
  static boost::mt19937 rng(seed);  //Random Number Generator
   
  boost::normal_distribution<double> gaussian_dist(mean, stdDev); //Normal Distribution
 
  // Create a Gaussian Random Number generator by binding 
  // with previously defined normal distribution object   
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > generator(rng, gaussian_dist);
 
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) 
  {
    // sample from the distribution
    (*this)[i] = generator();    
  }
  
  return;
};

// -------------------------------------------------
//TODO : find a smart way to define seed
void uqTeuchosVectorClass::cwSetGaussian(const uqTeuchosVectorClass& meanVec, const uqTeuchosVectorClass& stdDevVec)
{
  int seed =1;
  static boost::mt19937 rng(seed);  //Random Number Generator
     
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) 
  {
    boost::normal_distribution<double> gaussian_dist(meanVec[i], stdDevVec[i]); //Normal Distribution
 
  // Create a Gaussian Random Number generator by binding 
  // with previously defined normal distribution object   
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > generator(rng, gaussian_dist);
 
  // sample from the distribution
    (*this)[i] = generator();    
  }
  
  return;
};
// -------------------------------------------------
//TODO : find a smart way to define seed
void uqTeuchosVectorClass::cwSetGamma(const uqTeuchosVectorClass& aVec, const uqTeuchosVectorClass& bVec)
{
  int seed =1;
  static boost::mt19937 rng(seed);  //Random Number Generator
   
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) 
  {
    boost::gamma_distribution<double> gamma_dist(aVec[i], bVec[i]);   // Choose Gamma Distribution 
 
  // Create a Gamma Random Number generator by binding 
  // with previously defined normal distribution object   
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<double> > generator(rng, gamma_dist);
 
  // sample from the distribution
    (*this)[i] = generator();    
  }
  
  return;
};


// -------------------------------------------------
//TODO : find a smart way to define seed
// Using Gamma Distribution to calculate InverseGamma.
// Note the divisions: 1.0/b and the 1.0/generator; they are crucial
void uqTeuchosVectorClass::cwSetInverseGamma(const uqTeuchosVectorClass& aVec, const uqTeuchosVectorClass& bVec)
{
  int seed =1;
  static boost::mt19937 rng(seed);  //Random Number Generator
   
  for (unsigned int i = 0; i < this->sizeLocal(); ++i) 
  {
    // Choose Gamma Distribution with parabmer 1.0/b
    boost::gamma_distribution<double> gamma_dist(aVec[i], 1.0/bVec[i]);   
 
  // Create a Gamma Random Number generator by binding 
  // with previously defined normal distribution object   
  // sample from the distribution
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<double> > generator(rng, gamma_dist);
 
  // return 1.0 of the draw value
    (*this)[i] = 1.0/generator();    
  }
  
  return;
};

// -------------------------------------------------
//TODO : find a smart way to define seed
void uqTeuchosVectorClass::cwSetBeta(const uqTeuchosVectorClass& alpha, const uqTeuchosVectorClass& beta)
{
  int seed =1;  

  for (unsigned int i = 0; i < this->sizeLocal(); ++i) 
  {
   // draw a random number from (0,1)
   double randFromUnif = GetRandomDoubleUsingUniformZeroOneDistribution(seed); 
   // Choose Beta distribution with parameters alpha and beta
   boost::math::beta_distribution<double> beta_dist(alpha[i], beta[i]); 
   
   double randFromDist = quantile(beta_dist, randFromUnif);
   (*this)[i] = randFromDist;
  }
  
  return;
};

// -------------------------------------------------
void
uqTeuchosVectorClass::copy_to_std_vector(std::vector<double>& vec)
{
  unsigned int size = this->sizeLocal();
  vec.resize(size);
  
   for (unsigned int i = 0; i < size; ++i) 
  {
    vec[i] = m_vec[i];
  }
  
  return;
}

// -------------------------------------------------
// TODO: improve the message and instructions inside if.
void
uqTeuchosVectorClass::copy_from_std_vector(const std::vector<double> vec)
{
  
  unsigned int size = vec.size(), size2= this->sizeLocal();
  if (size != size2)
  {
    cout << endl << "Vectors have different size!" << endl;
  }
  
   for (unsigned int i = 0; i < size; ++i) 
  {
   m_vec[i] = vec[i];
  }
  
  return;
}

// -------------------------------------------------
void uqTeuchosVectorClass::sort()
{
  std::vector<double> vec;
  
  (*this).copy_to_std_vector(vec);  
  
  // using default comparison (operator <):
  std::sort (vec.begin(), vec.end()); 
  
  (*this).copy_from_std_vector(vec); 
};
/*====================================================*/

// operators ---------------------------------------
// element access ----------------------------------
double& uqTeuchosVectorClass::operator[](unsigned int i)
{
  return m_vec[i];
}


const double& uqTeuchosVectorClass::operator[](unsigned int i) const
{
  return m_vec[i];
}



//-------------------------------------------------
// assigns values a to all entrances in the vector
uqTeuchosVectorClass& uqTeuchosVectorClass::operator=(double a)
{
  m_vec.putScalar(a);
  
  return *this;
}



//-------------------------------------------------
uqTeuchosVectorClass& uqTeuchosVectorClass::operator=(const uqTeuchosVectorClass& rhs)
{
  unsigned int size1 = m_vec.length();
  unsigned int size2 = rhs.sizeLocal();
  
  UQ_FATAL_TEST_MACRO(size1 != size2, // mox
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::operator=()",
                      "the vectors do NOT have the same size.\n");
                      
  if (size1==size2){
    for (unsigned int i=0;i<size1;i++){
    m_vec[i]=rhs[i];  
    }
  }

  return *this;
}



//-------------------------------------------------
uqTeuchosVectorClass& uqTeuchosVectorClass::operator*=(double a)
{
  m_vec.scale(a); 	//Scale this vector by a; *this = a*this. 
                    
  return *this;
}

 //-------------------------------------------------
uqTeuchosVectorClass& uqTeuchosVectorClass::operator/=(double a)
{
  (*this) *= (1.0/a);

  return *this;
}

//-------------------------------------------------
uqTeuchosVectorClass& uqTeuchosVectorClass::operator*=(const uqTeuchosVectorClass& rhs)
{
  
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  UQ_FATAL_TEST_MACRO((size1 != size2),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::operator*=()",
                      "the vectors do NOT have the same size.\n");
                      
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] *= rhs[i];
    }
 }
 
  
  return *this;
}

//-------------------------------------------------
uqTeuchosVectorClass& uqTeuchosVectorClass::operator/=(const uqTeuchosVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  
    UQ_FATAL_TEST_MACRO((size1 != size2),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::operator/=()",
                      "the vectors do NOT have the same size.\n");
                      
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] /= rhs[i];
    }
  }
 
  
  return *this;
}

//-------------------------------------------------
uqTeuchosVectorClass& uqTeuchosVectorClass::operator+=(const uqTeuchosVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] += rhs[i];
    }
  }
  else{
    std::cout << "ERROR in uqTeuchosVectorClass:: operator+=()  ---> the vectors do NOT have the same size.\n";
  }
  
  return *this;
}



//-------------------------------------------------
uqTeuchosVectorClass& uqTeuchosVectorClass::operator-=(const uqTeuchosVectorClass& rhs)
{
  unsigned int size1 = this->sizeLocal();
  unsigned int size2 = rhs.sizeLocal();
  
  if (size1==size2){
    for (unsigned int i = 0; i < size1; ++i) {
      (*this)[i] -= rhs[i];
    }
  }
  else{
    std::cout << "ERROR in uqTeuchosVectorClass:: operator-=()  ---> the vectors do NOT have the same size.\n";
  }
  
  return *this;
}


//-------------------------------------------------
std::ostream&
operator<<(std::ostream& os, const uqTeuchosVectorClass& obj)
{
  obj.print(os);

  return os;
}

//-------------------------------------------------
unsigned int uqTeuchosVectorClass::sizeLocal() const
{
  UQ_FATAL_TEST_MACRO(m_vec.length() != (int) m_map.NumMyElements(),
                       m_env.worldRank(),
                       "uqTeuchosVectorClass::sizeLocal()",
                       "incompatible vec size");
  
  return m_vec.length();
}

//-------------------------------------------------
unsigned int uqTeuchosVectorClass::sizeGlobal() const
{
   UQ_FATAL_TEST_MACRO(m_vec.length() != (int) m_map.NumGlobalElements(),
                       m_env.worldRank(),
                       "uqTeuchosVectorClass::sizeGlobal()",
                       "incompatible vec size");
  return m_vec.length();
}

//-------------------------------------------------
double*
uqTeuchosVectorClass::values()
{
  return  m_vec.values();
};


//-------------------------------------------------
// Norms ------------------------------------------
double uqTeuchosVectorClass::norm2Sq() const
{  
 return (m_vec).dot(m_vec );  
}

double uqTeuchosVectorClass::norm2() const
{
  return std::sqrt(this->norm2Sq());
}


double uqTeuchosVectorClass::norm1() const
{
  double result = 0.;

  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += fabs((*this)[i]);
  }

  return result;
}


double uqTeuchosVectorClass::normInf() const
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
// ---------------------------------------------

double uqTeuchosVectorClass::sumOfComponents() const
{
  double result = 0.;
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    result += (*this)[i];
  }

  return result;
}


// Getting Max and Min values -------------------
// Max ------------------------------------------
double uqTeuchosVectorClass::getMaxValue( ) const  //dummy
{
  const unsigned int size = this->sizeLocal();
  std::vector<double> aux;
    
  for (unsigned int i=0; i<size; i++ ) {
    aux.push_back((*this)[i]) ;
  }
 
  return *max_element (aux.begin(),aux.end());
}

// Min ------------------------------------------
double uqTeuchosVectorClass::getMinValue( ) const //dummy
{  
   const unsigned int size = this->sizeLocal();
   std::vector<double> aux;
  
  for (unsigned int i=0; i<size; i++ ) {
    aux.push_back((*this)[i]) ;
    }
  
 return *min_element (aux.begin(),aux.end());
  
}

// Max index -----------------------------------
int uqTeuchosVectorClass::getMaxValueIndex( ) const   //dummy
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
int uqTeuchosVectorClass::getMinValueIndex( ) const //dummy
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
void uqTeuchosVectorClass::getMaxValueAndIndex( double& max_value, int& max_value_index ) //dummy
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
void uqTeuchosVectorClass::getMinValueAndIndex( double& min_value, int& min_value_index ) //dummy
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

// absolute values -------------------------------
// tested 1/8/13
uqTeuchosVectorClass
uqTeuchosVectorClass::abs() const
{
  uqTeuchosVectorClass abs_of_this_vec( *this );

  unsigned int size = abs_of_this_vec.sizeLocal();

  for( unsigned int i = 0; i < size; ++i )
    {
      abs_of_this_vec[i] = std::fabs( (*this)[i] );
    }

  return abs_of_this_vec;

}

//--------------------------------------------------------
void
uqTeuchosVectorClass::copy(const uqTeuchosVectorClass& rhs)
{
  this->uqVectorClass::copy(rhs); 
  
  unsigned int size1 = m_vec.length();
  unsigned int size2 = rhs.sizeLocal();
 
  if (size1==size2){
    for (unsigned int i=0;i<size1;i++){
    m_vec[i]=rhs[i];  
    }
  }
  return;
}

//====================================================
uqTeuchosVectorClass operator/(double a, const uqTeuchosVectorClass& x)
{
  uqTeuchosVectorClass answer(x); //copy x to answer
  
  answer.cwInvert();
  answer *= a;

  return answer;
}

uqTeuchosVectorClass operator/(const uqTeuchosVectorClass& x, const uqTeuchosVectorClass& y)
{
  uqTeuchosVectorClass answer(x);
  answer /= y;

  return answer;
}

uqTeuchosVectorClass operator*(double a, const uqTeuchosVectorClass& x)
{
  uqTeuchosVectorClass answer(x);
  answer *= a;

  return answer;
}

uqTeuchosVectorClass operator*(const uqTeuchosVectorClass& x, const uqTeuchosVectorClass& y)
{
  uqTeuchosVectorClass answer(x);
  answer *= y;

  return answer;
}

double scalarProduct(const uqTeuchosVectorClass& x, const uqTeuchosVectorClass& y)
{
  unsigned int size1 = x.sizeLocal();
  
//   UQ_FATAL_TEST_MACRO((size1 != size2),
//                       x.env().worldRank(),
//                       "scalarProduct()",
//                       "different sizes of x and y");

  
  double result = 0.;
  for (unsigned int i = 0; i < size1; ++i) {
    result += x[i]*y[i];
  }

   return result;
}

uqTeuchosVectorClass operator+(const uqTeuchosVectorClass& x, const uqTeuchosVectorClass& y)
{
  uqTeuchosVectorClass answer(x);
  answer += y;

  return answer;
}

uqTeuchosVectorClass operator-(const uqTeuchosVectorClass& x, const uqTeuchosVectorClass& y)
{
  uqTeuchosVectorClass answer(x);
  answer -= y;

  return answer;
}

bool
operator== (const uqTeuchosVectorClass& lhs, const uqTeuchosVectorClass& rhs)
{
  bool answer = true;

  unsigned int size1 = lhs.sizeLocal();
//  unsigned int size2 = rhs.sizeLocal();
//   UQ_FATAL_TEST_MACRO((size1 != size2),
//                       lhs.env().worldRank(),
//                       "operator==()",
//                       "different sizes of lhs and rhs");

  for (unsigned int i = 0; i < size1; ++i) {
    if (lhs[i] != rhs[i]) {
      answer = false;
      break;
    }
  }

  return answer;
}



//----------------------------------------------------
//----------------------------------------------------
void uqTeuchosVectorClass::cwSet(double value)
{
  (*this)=value;
  
  return;
}

//----------------------------------------------------
void uqTeuchosVectorClass::cwSet(unsigned int initialPos, const uqTeuchosVectorClass& vec)
{
//   UQ_FATAL_TEST_MACRO(initialPos >= this->sizeLocal(),
//                       m_env.worldRank(),
//                       "uqTeuchosVectorClass::cwSet()",
//                       "invalid initialPos");
// 
//   UQ_FATAL_TEST_MACRO((initialPos +vec.sizeLocal()) > this->sizeLocal(),
//                       m_env.worldRank(),
//                       "uqTeuchosVectorClass::cwSet()",
//                       "invalid vec.sizeLocal()");
  
  if (initialPos >= this->sizeLocal()){
    cout<< "ERROR in uqTeuchosVectorClass:: cwSet()  ---> invalid initialPos.\n";
    return;
  }

  if ((initialPos +vec.sizeLocal()) > this->sizeLocal()){
    cout<< "ERROR in uqTeuchosVectorClass:: cwSet()  ---> invalid invalid vec.sizeLocal().\n";
    return;
  }
  
  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    (*this)[initialPos+i] = vec[i];
  }

  return;
}


//----------------------------------------------------
/* extract elements from vector (*this) starting at position initialPos and save in vec */
void uqTeuchosVectorClass::cwExtract(unsigned int initialPos, uqTeuchosVectorClass& vec) const
{
//   UQ_FATAL_TEST_MACRO(initialPos >= this->sizeLocal(),
//                       m_env.worldRank(),
//                       "uqTeuchosVectorClass::cwExtract()",
//                       "invalid initialPos");
// 
//   UQ_FATAL_TEST_MACRO((initialPos +vec.sizeLocal()) > this->sizeLocal(),
//                       m_env.worldRank(),
//                       "uqTeuchosVectorClass::cwExtract()",
//                       "invalid vec.sizeLocal()");

 if (initialPos >= this->sizeLocal()){
    cout<< "ERROR in uqTeuchosVectorClass:: cwExtract()  ---> invalid initialPos.\n";
    return;
  }

  if ((initialPos +vec.sizeLocal()) > this->sizeLocal()){
    cout<< "ERROR in uqTeuchosVectorClass:: cwExtract()  ---> invalid invalid vec.sizeLocal().\n";
    return;
  }
  
  for (unsigned int i = 0; i < vec.sizeLocal(); ++i) {
    vec[i] = (*this)[initialPos+i];
  }

  return;
}


void uqTeuchosVectorClass::cwInvert()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = 1./(*this)[i];
  }

  return;
}

void uqTeuchosVectorClass::cwSqrt()
{
  unsigned int size = this->sizeLocal();
  for (unsigned int i = 0; i < size; ++i) {
    (*this)[i] = sqrt((*this)[i]);
  }

  return;
}


/*----------------------------------------------------------------
*  DLARNV returns a vector of n random real numbers from a uniform or
*  normal distribution.
*
*  Arguments
*  =========
*  IDIST   (input) INTEGER
*          Specifies the distribution of the random numbers:
*          = 1:  uniform (0,1)
*          = 2:  uniform (-1,1)
*          = 3:  normal (0,1)
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be odd.
*          On exit, the seed is updated.
*  N       (input) INTEGER
*          The number of random numbers to be generated.
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          The generated random numbers.
* Also view: 
* http://www.cs.cmu.edu/~dlr/dlr_libs/dlrlibs1.28/html/pseudoRandom_8cpp-source.html
*------------------------------------------------------------------------*/
void uqTeuchosVectorClass::cwSetGaussian2(double mean, double stdDev)
{
  unsigned int iseed_size = 4; // Lapack expects iseed_size = 4, always.
  int iseed[iseed_size];
  Teuchos::SerialDenseVector<int,double> x(iseed_size);
  
  x.random();
  x.scale(10);
    
  for(unsigned int i=0; i<iseed_size; i++)
  {
    iseed[i] = fabs(floor(x[i])); 
  }
   // iseed[3] must be odd.
  if(iseed[3] % 2 == 0) {
    iseed[3] -= 1;
  }

  Teuchos::LAPACK<int, double> lapack;
  int idist=3;  	     // Tells lapack we want a normal distribution
  int dim=this->sizeLocal(); // Amount of values needed
  double returnValue[dim];   // The number we'll return
     
  lapack.LARNV(idist,iseed,dim,returnValue); 	
  
  for (int i=0; i < dim; i++){
         (*this)[i] = mean + (returnValue[i] * stdDev);
   }
   
  return;
}

//---------------------------------------------------- 
 void uqTeuchosVectorClass::cwSetGaussian2(const uqTeuchosVectorClass& meanVec, const uqTeuchosVectorClass& stdDevVec)
{
  unsigned int iseed_size = 4; // Lapack expects iseed_size = 4, always.
  int iseed[iseed_size];
  Teuchos::SerialDenseVector<int,double> x(iseed_size);
  
  x.random();
  x.scale(10);
    
  for(unsigned int i=0; i<iseed_size; i++)
  {
    iseed[i] = fabs(floor(x[i])); 
  }
   // iseed[3] must be odd.
  if(iseed[3] % 2 == 0) {
    iseed[3] -= 1;
  }

  Teuchos::LAPACK<int, double> lapack;
  int idist=3;  	     // Tells lapack we want a normal distribution
  int dim=this->sizeLocal(); // Amount of values needed
  double returnValue[dim];   // The number we'll return
     
  lapack.LARNV(idist,iseed,dim,returnValue); 	
  
  for (int i=0; i < dim; i++){
         (*this)[i] = meanVec[i] + (returnValue[i] * stdDevVec[i]);
   }
   
  return;
}


//----------------------------------------------------
 void uqTeuchosVectorClass::cwSetUniform2(const uqTeuchosVectorClass& lowerBoundaVec, const uqTeuchosVectorClass& upperBoundVec)
{
  unsigned int iseed_size = 4; // Lapack expects iseed_size = 4, always.
  int iseed[iseed_size];
  Teuchos::SerialDenseVector<int,double> x(iseed_size);
  
  x.random();
  x.scale(10);
    
  for(unsigned int i=0; i<iseed_size; i++)
  {
    iseed[i] = fabs(floor(x[i])); 
  }
   // iseed[3] must be odd.
  if(iseed[3] % 2 == 0) {
    iseed[3] -= 1;
  }

  Teuchos::LAPACK<int, double> lapack;
  int idist=1;  	     // Tells lapack we want a uniform distribution
  int dim=this->sizeLocal(); // Amount of values needed
  double returnValue[dim];   // The number we'll return
     
  lapack.LARNV(idist,iseed,dim,returnValue); 	
  
  for (int i=0; i < dim; i++){
         (*this)[i] = lowerBoundaVec[i] + (upperBoundVec[i]-lowerBoundaVec[i])*returnValue[i];
   }
   
  return;
}



//----------------------------------------------------
//Boost version
// implemented/checked 2/26/13
//TODO : find better way to get the seed
 void uqTeuchosVectorClass::cwSetUniform(const uqTeuchosVectorClass& lowerBoundaVec, const uqTeuchosVectorClass& upperBoundVec)
{
   int seed = 1;

  for (unsigned int i=0; i < this->sizeLocal(); i++){
    double random= GetRandomDoubleUsingUniformZeroOneDistribution(seed);
    (*this)[i] = lowerBoundaVec[i] + (upperBoundVec[i]-lowerBoundaVec[i])*random;
   }
  return;
}
//----------------------------------------------------

void
uqTeuchosVectorClass::cwSetConcatenated(const uqTeuchosVectorClass& v1, const uqTeuchosVectorClass& v2)
{
//   UQ_FATAL_TEST_MACRO(this->sizeLocal() != v1.sizeLocal() + v2.sizeLocal(),
//                       m_env.worldRank(),
//                       "uqTeuchosVectorClass::cwSetConcatenated(1)",
//                       "incompatible vector sizes");

  if ( this->sizeLocal() != v1.sizeLocal() + v2.sizeLocal() ) {
    std::cout << "ERROR in uqTeuchosVectorClass:: cwSetConcatenated  ---> the vectors' sizes are not compatible.\n";
    cout<< "      in uqTeuchosVectorClass:: cwSetConcatenated  ---> resizing resulting vector... new size = " << v1.sizeLocal()+v1.sizeLocal() <<endl;
    
    m_vec.resize(v1.sizeLocal()+v1.sizeLocal());
  }
   
  for (unsigned int i = 0; i < v1.sizeLocal(); ++i) {
    (*this)[i] = v1[i];
  }

  for (unsigned int i = 0; i < v2.sizeLocal(); ++i) {
    (*this)[v1.sizeLocal()+i] = v2[i];
  }

  return;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++

bool
uqTeuchosVectorClass::atLeastOneComponentSmallerThan(const uqTeuchosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::atLeastOneComponentSmallerThan()",
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
uqTeuchosVectorClass::atLeastOneComponentBiggerThan(const uqTeuchosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::atLeastOneComponentBiggerThan()",
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

bool
uqTeuchosVectorClass::atLeastOneComponentSmallerOrEqualThan(const uqTeuchosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::atLeastOneComponentSmallerOrEqualThan()",
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
uqTeuchosVectorClass::atLeastOneComponentBiggerOrEqualThan(const uqTeuchosVectorClass& rhs) const
{
  UQ_FATAL_TEST_MACRO((this->sizeLocal() != rhs.sizeLocal()),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::atLeastOneComponentBiggerOrEqualThan()",
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

//++++++++++++++++++++++++++++++++++++++++++++++

void
uqTeuchosVectorClass::print(std::ostream& os) const
{
  //std::cout << "In uqTeuchosVectorClass::print(): before sizelocal()"
  //          << std::endl;
  unsigned int size = this->sizeLocal();
  //std::cout << "In uqTeuchosVectorClass::print(): after sizelocal()"
  //          << std::endl;

  //std::cout << "In uqTeuchosVectorClass::print(): before os.flags()"
  //          << std::endl;
  std::ostream::fmtflags curr_fmt = os.flags();
  //std::cout << "In uqTeuchosVectorClass::print(): after os.flags()"
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
      //std::cout << "In uqTeuchosVectorClass::print(): where expected"
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

  //std::cout << "In uqTeuchosVectorClass::print(): before os.flags(curr_fmt)"
  //          << std::endl;
  os.flags(curr_fmt);
  //std::cout << "In uqTeuchosVectorClass::print(): after os.flags(curr_fmt)"
  //          << std::endl;

  return;
}
//------------------------------------------------------------

void
uqTeuchosVectorClass::subReadContents(
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds)
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::subReadContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::subReadContents()",
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
                        "uqTeuchosVectorClass::subReadContents()",
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
                          "uqTeuchosVectorClass::subReadContents()",
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
                          "uqTeuchosVectorClass::subReadContents()",
                          "symbol ')' not found in first line of file");
      nParamsString[posInParamsString++] = tmpString[posInTmpString++];
    } while (tmpString[posInTmpString] != ')');
    nParamsString[posInParamsString] = '\0';

    // Convert 'n_positions' and 'n_params' strings to numbers
    unsigned int sizeOfVecInFile = (unsigned int) strtod(nPositionsString,NULL);
    unsigned int numParamsInFile = (unsigned int) strtod(nParamsString,   NULL);
    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqTeuchosVectorClass::subReadContents()"
                              << ": fullRank "            << m_env.fullRank()
                              << ", sizeOfVecInFile = "   << sizeOfVecInFile
                              << ", numParamsInFile = "   << numParamsInFile
                              << ", this->sizeLocal() = " << this->sizeLocal()
                              << std::endl;
    }

    // Check if [size of vec in file] >= [requested sub vec size]
    UQ_FATAL_TEST_MACRO(sizeOfVecInFile < subReadSize,
                        m_env.worldRank(),
                        "uqTeuchosVectorClass::subReadContents()",
                        "size of vec in file is not big enough");

    // Check if [num params in file] == [num params in current vec]
    UQ_FATAL_TEST_MACRO(numParamsInFile != numParams,
                        m_env.worldRank(),
                        "uqTeuchosVectorClass::subReadContents()",
                        "number of parameters of vec in file is different than number of parameters in this vec object");

    // Code common to any core in a communicator
    unsigned int maxCharsPerLine = 64*numParams; // Up to about 60 characters to represent each parameter value

    unsigned int lineId = 0;
    while (lineId < idOfMyFirstLine) {
      filePtrSet.ifsVar->ignore(maxCharsPerLine,'\n');
      lineId++;
    };

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqTeuchosVectorClass::subReadContents()"
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
                        "uqTeuchosVectorClass::subReadContents()",
                        "in core 0, string should be the '=' sign");

    // Take into account the ' [' portion
    std::streampos tmpPos = filePtrSet.ifsVar->tellg();
    filePtrSet.ifsVar->seekg(tmpPos+(std::streampos)2);

    if (m_env.subDisplayFile()) {
      *m_env.subDisplayFile() << "In uqTeuchosVectorClass::subReadContents()"
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

//------------------------------------------------------------------

void
uqTeuchosVectorClass::subWriteContents(
  const std::string&            varNamePrefix,
  const std::string&            fileName,
  const std::string&            fileType,
  const std::set<unsigned int>& allowedSubEnvIds) const
{
  UQ_FATAL_TEST_MACRO(m_env.subRank() < 0,
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::subWriteContents()",
                      "unexpected subRank");

  UQ_FATAL_TEST_MACRO(this->numOfProcsForStorage() > 1,
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::subWriteContents()",
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
//------------------------------------------------------------------
// added 2/28/13
void
uqTeuchosVectorClass::mpiBcast(int srcRank, const uqMpiCommClass& bcastComm)
{
  // Filter out those nodes that should not participate
  if (bcastComm.MyPID() < 0) return;

  // Check 'srcRank'
  UQ_FATAL_TEST_MACRO((srcRank < 0) || (srcRank >= bcastComm.NumProc()),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::mpiBcast()",
                      "invalud srcRank");

  // Check number of participant nodes
  double localNumNodes = 1.;
  double totalNumNodes = 0.;
  bcastComm.Allreduce((void *) &localNumNodes, (void *) &totalNumNodes, (int) 1, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                      "uqTeuchosVectorClass::mpiBcast()",
                      "failed MPI.Allreduce() for numNodes");
  UQ_FATAL_TEST_MACRO(((int) totalNumNodes) != bcastComm.NumProc(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::mpiBcast()",
                      "inconsistent numNodes");

  // Check that all participant nodes have the same vector size
  double localVectorSize  = this->sizeLocal();
  double sumOfVectorSizes = 0.; 
  bcastComm.Allreduce((void *) &localVectorSize, (void *) &sumOfVectorSizes, (int) 1, uqRawValue_MPI_DOUBLE, uqRawValue_MPI_SUM,
                      "uqTeuchosVectorClass::mpiBcast()",
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
                      "uqTeuchosVectorClass::mpiBcast()",
                      "inconsistent vectorSize");

  // Ok, bcast data
  std::vector<double> dataBuffer((unsigned int) localVectorSize, 0.);
  if (bcastComm.MyPID() == srcRank) {
    for (unsigned int i = 0; i < dataBuffer.size(); ++i) {
      dataBuffer[i] = (*this)[i];
    }
  }

  bcastComm.Bcast((void *) &dataBuffer[0], (int) localVectorSize, uqRawValue_MPI_DOUBLE, srcRank,
                  "uqTeuchosVectorClass::mpiBcast()",
                  "failed MPI.Bcast()");

  if (bcastComm.MyPID() != srcRank) {
    for (unsigned int i = 0; i < dataBuffer.size(); ++i) {
      (*this)[i] = dataBuffer[i];
    }
  }

  return;
}
//------------------------------------------------------------------
// added 2/28/13

void
uqTeuchosVectorClass::mpiAllReduce(uqRawType_MPI_Op mpiOperation, const uqMpiCommClass& opComm, uqTeuchosVectorClass& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  unsigned int size = this->sizeLocal();
  UQ_FATAL_TEST_MACRO(size != resultVec.sizeLocal(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::mpiAllReduce()",
                      "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double srcValue = (*this)[i];
    double resultValue = 0.;
    opComm.Allreduce((void *) &srcValue, (void *) &resultValue, (int) 1, uqRawValue_MPI_DOUBLE, mpiOperation,
                     "uqTeuchosVectorClass::mpiAllReduce()",
                     "failed MPI.Allreduce()");
    resultVec[i] = resultValue;
  }

  return;
}
//------------------------------------------------------------------
// added 2/28/13

void
uqTeuchosVectorClass::mpiAllQuantile(double probability, const uqMpiCommClass& opComm, uqTeuchosVectorClass& resultVec) const
{
  // Filter out those nodes that should not participate
  if (opComm.MyPID() < 0) return;

  UQ_FATAL_TEST_MACRO((probability < 0.) || (1. < probability),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::mpiAllQuantile()",
                      "invalid input");

  unsigned int size = this->sizeLocal();
  UQ_FATAL_TEST_MACRO(size != resultVec.sizeLocal(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::mpiAllQuantile()",
                      "different vector sizes");

  for (unsigned int i = 0; i < size; ++i) {
    double auxDouble = (int) (*this)[i];
    std::vector<double> vecOfDoubles(opComm.NumProc(),0.);
    opComm.Gather((void *) &auxDouble, 1, uqRawValue_MPI_DOUBLE, (void *) &vecOfDoubles[0], (int) 1, uqRawValue_MPI_DOUBLE, 0,
                  "uqTeuchosVectorClass::mpiAllQuantile()",
                  "failed MPI.Gather()");

    std::sort(vecOfDoubles.begin(), vecOfDoubles.end());

    double result = vecOfDoubles[(unsigned int)( probability*((double)(vecOfDoubles.size()-1)) )];

    opComm.Bcast((void *) &result, (int) 1, uqRawValue_MPI_DOUBLE, 0,
                 "uqTeuchosVectorClass::mpiAllQuantile()",
                 "failed MPI.Bcast()");

    resultVec[i] = result;
  }

  return;
}



//------------------------------------------------------------------
// added/tested 2/28/13
void
uqTeuchosVectorClass::matlabLinearInterpExtrap(
  const uqTeuchosVectorClass& x1Vec,
  const uqTeuchosVectorClass& y1Vec,
  const uqTeuchosVectorClass& x2Vec)
{
  UQ_FATAL_TEST_MACRO(x1Vec.sizeLocal() <= 1,
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::matlabLinearInterpExtrap()",
                      "invalid 'x1' size");

  UQ_FATAL_TEST_MACRO(x1Vec.sizeLocal() != y1Vec.sizeLocal(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::matlabLinearInterpExtrap()",
                      "invalid 'x1' and 'y1' sizes");

  UQ_FATAL_TEST_MACRO(x2Vec.sizeLocal() != this->sizeLocal(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::matlabLinearInterpExtrap()",
                      "invalid 'x2' and 'this' sizes");

  for (unsigned int i = 1; i < x1Vec.sizeLocal(); ++i) { // Yes, '1'
    UQ_FATAL_TEST_MACRO(x1Vec[i] <= x1Vec[i-1],
                        m_env.worldRank(),
                        "uqTeuchosVectorClass::matlabLinearInterpExtrap()",
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

//------------------------------------------------------------------
// added 2/28/13
void
uqTeuchosVectorClass::matlabDiff(
  unsigned int      firstPositionToStoreDiff,
  double            valueForRemainderPosition,
  uqTeuchosVectorClass& outputVec) const
{
  unsigned int size = this->sizeLocal();

  UQ_FATAL_TEST_MACRO(firstPositionToStoreDiff > 1,
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::matlabDiff()",
                      "invalid firstPositionToStoreDiff");

  UQ_FATAL_TEST_MACRO(size != outputVec.sizeLocal(),
                      m_env.worldRank(),
                      "uqTeuchosVectorClass::matlabDiff()",
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


//------------------------------------------------------------------
//------------------------------------------------------------------
double 
uqTeuchosVectorClass::GetRandomDoubleUsingNormalDistribution(int seed, double mean,double sigma)
{ 
  /* Choose RNG: mt19937 is the second fastest (93%) and 
   * generates good uniform distribution in up to 623 dimensions  */
  static boost::mt19937 rng(seed); 
   
  /* Choose Normal Distribution */
  boost::normal_distribution<double> gaussian_dist(mean, sigma);
 
  /* Create a Gaussian Random Number generator by binding 
   * with previously defined normal distribution object   */
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > generator(rng, gaussian_dist);
 
  // sample from the distribution
  return generator();
  
}
//------------------------------------------------------------------
double uqTeuchosVectorClass::GetRandomDoubleUsingUniformZeroOneDistribution(int seed)
{
  boost::mt19937 rng(seed);
  static boost::uniform_01<boost::mt19937> zeroone(rng);
  return zeroone();
}
#endif // ifdef QUESO_HAS_TRILINOS

