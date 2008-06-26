#ifndef __UQ_GSL_VECTOR_H__
#define __UQ_GSL_VECTOR_H__

#include <uqVector.h>
#include <gsl/gsl_vector.h>

class uqGslVectorClass : public uqVectorClass
{
public:
  uqGslVectorClass();
  uqGslVectorClass(const uqEnvironmentClass& env, unsigned int size);
  uqGslVectorClass(const uqEnvironmentClass& env, double d1, double d2, unsigned int size); // MATLAB linspace
  uqGslVectorClass(const uqGslVectorClass& y);
 ~uqGslVectorClass();

  uqGslVectorClass& operator= (const uqGslVectorClass& rhs);
  uqGslVectorClass& operator*=(double a);
  uqGslVectorClass& operator/=(double a);
  uqGslVectorClass& operator*=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator/=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator+=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator-=(const uqGslVectorClass& rhs);
            double& operator[](unsigned int i);
      const double& operator[](unsigned int i) const;

  unsigned int size             () const;
  double       norm2Sq          () const;
  double       norm2            () const;
  double       sumOfComponents  () const;
  void         cwSet            (double value);
  void         cwSetGaussian    (gsl_rng* rng, double mean, double stdDev);
  void         cwInvert         ();
  void         sort             ();
  void         print            (std::ostream& os) const;

  bool         atLeastOneComponentSmallerThan(const uqGslVectorClass& rhs) const;
  bool         atLeastOneComponentBiggerThan (const uqGslVectorClass& rhs) const;
  gsl_vector*  data             () const; // Necessary for uqGslMatrixClass::invertMultiply()

private:

  void         copy             (const uqGslVectorClass& src);

  gsl_vector* m_vec;
};

uqGslVectorClass operator/    (const double a,            const uqGslVectorClass& x  );
uqGslVectorClass operator/    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator*    (const double a,            const uqGslVectorClass& x  );
uqGslVectorClass operator*    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
double           scalarProduct(const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator+    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator-    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
std::ostream&    operator<<   (std::ostream& os,          const uqGslVectorClass& obj);

#endif // __UQ_GSL_VECTOR_H__
