#ifndef __UQ_MATRIX_H__
#define __UQ_MATRIX_H__

#include <uqEnvironment.h>
#include <uqVector.h>
#include <iostream>

class uqMatrixClass
{
public:
           uqMatrixClass();
           uqMatrixClass(const uqEnvironmentClass& env);
           uqMatrixClass(const uqMatrixClass& rhs);
  virtual ~uqMatrixClass();

  uqMatrixClass& operator= (const uqMatrixClass& rhs);
  uqMatrixClass& operator*=(double a);
  uqMatrixClass& operator+=(const uqMatrixClass& rhs);
  uqMatrixClass& operator-=(const uqMatrixClass& rhs);

  const uqEnvironmentClass& env() const;

  virtual unsigned int   numRows       () const = 0;
  virtual unsigned int   numCols       () const = 0;
  virtual int            chol          () = 0;
  virtual void           zeroLower     (bool includeDiagonal = false) = 0;
  virtual void           zeroUpper     (bool includeDiagonal = false) = 0;
  virtual void           print         (std::ostream& os) const = 0;

protected:
  virtual void           copy          (const uqMatrixClass& src);

  const uqEnvironmentClass& m_env;
};

#endif // __UQ_MATRIX_H__
