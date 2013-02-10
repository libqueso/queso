#ifndef __UQ_TEUCHOS_VECTOR_H__
#define __UQ_TEUCHOS_VECTOR_H__


#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_LAPACK.hpp>
#include <vector>
#include <uqVector.h>
#include <gsl/gsl_vector.h>
#include "boost/random.hpp"

class uqTeuchosVectorClass : public uqVectorClass
{
public:
  uqTeuchosVectorClass();
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map);
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, const uqMapClass& map, double value);
  uqTeuchosVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const uqMapClass& map); // MATLAB linspace
  uqTeuchosVectorClass(const uqTeuchosVectorClass&         v, double d1, double d2);                        // MATLAB linspace
  uqTeuchosVectorClass(const uqTeuchosVectorClass&         y);
  
  ~uqTeuchosVectorClass();

  uqTeuchosVectorClass abs() const;

  uqTeuchosVectorClass& operator= (double a);
  uqTeuchosVectorClass& operator= (const uqTeuchosVectorClass& rhs);
  uqTeuchosVectorClass& operator*=(double a);
  uqTeuchosVectorClass& operator/=(double a);
  uqTeuchosVectorClass& operator*=(const uqTeuchosVectorClass& rhs);
  uqTeuchosVectorClass& operator/=(const uqTeuchosVectorClass& rhs);
  uqTeuchosVectorClass& operator+=(const uqTeuchosVectorClass& rhs);
  uqTeuchosVectorClass& operator-=(const uqTeuchosVectorClass& rhs);
  

  double& operator[](unsigned int i);
  const double& operator[](unsigned int i) const;
  
  unsigned int sizeLocal        () const;
  unsigned int sizeGlobal       () const;
  double*      values           () ;//added by Kemelli on 12/05/12, precisa checar
  double       norm2Sq          () const;
  double       norm2            () const;
  double       norm1            () const;
  double       normInf          () const;
  double       sumOfComponents  () const;
  
  double       getMaxValue      () const;
  double       getMinValue      () const;
  int          getMaxValueIndex () const;
  int          getMinValueIndex () const;
  void         getMaxValueAndIndex( double& value, int& index );
  void         getMinValueAndIndex( double& value, int& index );
  
  bool         atLeastOneComponentSmallerThan       (const uqTeuchosVectorClass& rhs) const;
  bool         atLeastOneComponentBiggerThan        (const uqTeuchosVectorClass& rhs) const;
  bool         atLeastOneComponentSmallerOrEqualThan(const uqTeuchosVectorClass& rhs) const;
  bool         atLeastOneComponentBiggerOrEqualThan (const uqTeuchosVectorClass& rhs) const;
  
  void	       cwSet(double value);
  void         cwSet(unsigned int initialPos, const uqTeuchosVectorClass& vec);
  void         cwExtract(unsigned int initialPos, uqTeuchosVectorClass& vec) const; 
  void         cwInvert();
  void         cwSqrt();
  void	       cwSetGaussian(double mean, double stdDev);
  void	       cwSetGaussian(const uqTeuchosVectorClass& meanVec,       const uqTeuchosVectorClass& stdDevVec);
  void         cwSetUniform (const uqTeuchosVectorClass& lowerBoundVec, const uqTeuchosVectorClass& upperBoundVec);
  void         cwSetUniform2(const uqTeuchosVectorClass& lowerBoundVec, const uqTeuchosVectorClass& upperBoundVec);//dummy
   
  void         cwSetConcatenated(const uqTeuchosVectorClass& v1, const uqTeuchosVectorClass& v2);
  /*Nope*/
  /*
  void         cwSetBeta        (const uqTeuchosVectorClass& alpha,   const uqTeuchosVectorClass& beta     );
  void         cwSetInverseGamma(const uqTeuchosVectorClass& alpha,   const uqTeuchosVectorClass& beta     );
  */
  void         cwSetGaussian2   (double mean, double stdDev);
  void         cwSetGaussian2   (const uqTeuchosVectorClass& meanVec, const uqTeuchosVectorClass& stdDevVec);
  void         cwSetGamma       (const uqTeuchosVectorClass& a,       const uqTeuchosVectorClass& b        );
  void         sort             ();
  void         print            (std::ostream& os) const;
  void         subReadContents  (const std::string&            fileName,
				 const std::string&            fileType,
				 const std::set<unsigned int>& allowedSubEnvIds);
  void         subWriteContents (const std::string&            varNamePrefix,
                                 const std::string&            fileName,
                                 const std::string&            fileType,
                                 const std::set<unsigned int>& allowedSubEnvIds) const;
				 
  double	GetRandomDoubleUsingNormalDistribution(int seed, double mean,double sigma);
  void copy_to_std_vector(std::vector<double>& vec); // output
  void copy_from_std_vector(const std::vector<double> vec); //(input, output)


private:
  Teuchos::SerialDenseVector<int,double> m_vec;
  void         copy             (const uqTeuchosVectorClass& src);

};

uqTeuchosVectorClass copy	  (int , 			    double []);
uqTeuchosVectorClass operator/    (double a,            	    const uqTeuchosVectorClass& x  );
uqTeuchosVectorClass operator/    (const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
uqTeuchosVectorClass operator*    (double a,		            const uqTeuchosVectorClass& x  );
uqTeuchosVectorClass operator*    (const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
double               scalarProduct(const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
uqTeuchosVectorClass operator+    (const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
uqTeuchosVectorClass operator-    (const uqTeuchosVectorClass& x,   const uqTeuchosVectorClass& y  );
bool                 operator==   (const uqTeuchosVectorClass& lhs, const uqTeuchosVectorClass& rhs);
std::ostream&        operator<<   (std::ostream& os,                const uqTeuchosVectorClass& obj);

#endif //__UQ_TEUCHOS_VECTOR_H__


