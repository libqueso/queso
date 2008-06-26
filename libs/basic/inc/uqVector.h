#ifndef __UQ_VECTOR_H__
#define __UQ_VECTOR_H__

#include <uqEnvironment.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>

class uqVectorClass
{
public:
           uqVectorClass();
           uqVectorClass(const uqEnvironmentClass& env);
           uqVectorClass(const uqVectorClass& rhs);
  virtual ~uqVectorClass();

  uqVectorClass& operator= (const uqVectorClass& rhs);
  uqVectorClass& operator*=(double a);
  uqVectorClass& operator/=(double a);
  uqVectorClass& operator+=(const uqVectorClass& rhs);
  uqVectorClass& operator-=(const uqVectorClass& rhs);

    const uqEnvironmentClass& env                 ()           const;
          void                setPrintHorizontally(bool value) const; // Yes, 'const'
          bool                getPrintHorizontally()           const;

  virtual unsigned int        size                () const = 0;
  virtual void                cwSet               (double value) = 0;
  virtual void                cwSetGaussian       (gsl_rng* rng, double mean, double stdDev) = 0;
  virtual void                cwInvert            () = 0;
  virtual void                sort                () = 0;
  virtual void                print               (std::ostream& os) const = 0;

protected:
  virtual void                copy                (const uqVectorClass& src);

  const uqEnvironmentClass& m_env;
  mutable bool m_printHorizontally;
};

#endif // __UQ_VECTOR_H__
