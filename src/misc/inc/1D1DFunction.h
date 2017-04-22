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

#ifndef UQ_1D_1D_FUNCTION_H
#define UQ_1D_1D_FUNCTION_H

#include <queso/ScalarSequence.h>
#include <queso/1DQuadrature.h>
#include <queso/Environment.h>
#include <queso/Defines.h>
#include <vector>
#include <math.h>
#include <fstream>

namespace QUESO {

/*! \file 1D1DFunction.h
    \brief One-dimension function class.
*/



//*****************************************************
// Base 1D->1D class
//*****************************************************

/*! \class Base1D1DFunction
    \brief Class for one-dimensional functions.

    Base class for one-dimensional functions.
*/
class
Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  Base1D1DFunction(double minDomainValue,
			  double maxDomainValue);

  //! Destructor.
  virtual ~Base1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the minimum value of the domain of the (one-dimensional) function.
  double minDomainValue() const;

  //! Returns the maximum value of the domain of the (one-dimensional) function.
  double maxDomainValue() const;

  //! Returns the value of the (one-dimensional) function. See template specialization.
  virtual  double value         (double domainValue) const = 0;

  //! Returns the value of the derivative of the function. See template specialization.
  virtual  double deriv         (double domainValue) const = 0;

  //! TODO: Multiplies \c this function with \c function, and integrates it numerically.  See template specialization.
  /*! \todo: Please, implement me!*/
  virtual  double multiplyAndIntegrate(const Base1D1DFunction& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const;
  //@}
protected:
  double m_minDomainValue;
  double m_maxDomainValue;
};

//*****************************************************
// Generic 1D->1D class
//*****************************************************
/*! \class Generic1D1DFunction
    \brief Class for generic one-dimensional functions.

    This class implements generic one-dimensional functions.
*/
class Generic1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  Generic1D1DFunction(double minDomainValue,
                             double maxDomainValue,
                             double (*valueRoutinePtr)(double domainValue, const void* routinesDataPtr),
                             double (*derivRoutinePtr)(double domainValue, const void* routinesDataPtr),
                             const void* routinesDataPtr);
  //! Destructor.
  ~Generic1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the (one-dimensional) function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it evaluates the function at such point. */
  double value(double domainValue) const;

  //! Returns the value of the derivative of the function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it returns the value of the derivative at such point. */
  double deriv(double domainValue) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  double (*m_valueRoutinePtr)(double domainValue, const void* routinesDataPtr);
  double (*m_derivRoutinePtr)(double domainValue, const void* routinesDataPtr);
  const void* m_routinesDataPtr;
};

//*****************************************************
// Constant 1D->1D class
//*****************************************************
/*! \class Constant1D1DFunction
    \brief Class for constant one-dimensional functions.

    This class implements constant one-dimensional functions.
*/
class Constant1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  Constant1D1DFunction(double minDomainValue,
                              double maxDomainValue,
                              double constantValue);
  //! Destructor.
  ~Constant1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the constant function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it evaluates the function at such point, which is the constant
   * value \c constantValue passed to the constructor of this class. */
  double value(double domainValue) const;

  //! Returns the value of the derivative of the constant function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it returns 0 (derivative of a constant function. */
  double deriv(double domainValue) const;
//@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  double m_constantValue;
};

//*****************************************************
// Linear 1D->1D class
//*****************************************************
/*! \class Linear1D1DFunction
    \brief Class for linear one-dimensional functions.

    This class implements linear one-dimensional functions.
    A common linear function is the the equation of a line: \f[f(x) = y_1 + m (x - x_1), \f]
    which is a linear function with slope (or rate) \f$ m \f$ and passing through the point
    \f$(x_1,y_1)\f$.

    This class implements a linear function such as the one describe above where the rate is
    given by \c rateValue, and the point is <c>(referenceDomainValue,referenceImageValue)</c>.
*/
class Linear1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Constructs a linear function of constant rate of change \c rateValue which passes through
   *  the point (x,y)=<c>(referenceDomainValue,referenceImageValue)</c>.*/
  Linear1D1DFunction(double minDomainValue,
                            double maxDomainValue,
                            double referenceDomainValue,
                            double referenceImageValue,
                            double rateValue);
  //! Destructor.
  ~Linear1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the linear function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it evaluates the function at such point. */
  double value(double domainValue) const;

  //! Returns the value of the derivative of the linear function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it returns the value of the derivative at such point (namely
   * \c m_rateValue). */
  double deriv(double domainValue) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  //! Reference value in the function domain; \f$ x_1 \f$ in  \f$ f(x) = y_1 + m (x - x_1)\f$.
  double m_referenceDomainValue;

  //! Reference value in the function image; \f$ y_1 \f$ in  \f$ f(x) = y_1 + m (x - x_1)\f$.
  double m_referenceImageValue;

  //! Rate value; \f$ m\f$ in  \f$f(x) = y_1 + m (x - x_1)\f$.
  double m_rateValue;
};

//*****************************************************
// PiecewiseLinear 1D->1D class
//*****************************************************
/*! \class PiecewiseLinear1D1DFunction
    \brief Class for piecewise-linear one-dimensional functions.

    This class implements piecewise-linear one-dimensional functions.
*/
class PiecewiseLinear1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  PiecewiseLinear1D1DFunction(double                     minDomainValue,
                                     double                     maxDomainValue,
                                     const std::vector<double>& referenceDomainValues,
                                     double                     referenceImageValue0,
                                     const std::vector<double>& rateValues);
  //! Destructor.
  ~PiecewiseLinear1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the piecewise-linear function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it evaluates the function at such point. */
  double value(double domainValue) const;

  //! Returns the value of the derivative of the piecewise-linear function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it returns the value of the derivative at such point. */
  double deriv(double domainValue) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  //! Number of points which will be evaluated.
  unsigned int        m_numRefValues;

  //! Reference values in the piecewise-linear function domain; \f$ x_1i \f$ in   \f$ f_i(x) = y_{1i} + m_i (x - x_{1i})\f$, for each \f$ i \f$=1 .. \c m_numRefValues.
  std::vector<double> m_referenceDomainValues;

  //! Reference values in the piecewise-linear function image; \f$ y_{1i} \f$ in  \f$ f_i(x) = y_{1i} + m_i (x - x_{1i})\f$, for each \f$ i \f$=1 .. \c m_numRefValues.
  std::vector<double> m_referenceImageValues;

  //! Rate value; \f$ m_i \f$ in \f$ f_i(x) = y_{1i} + m_i (x - x_{1i})\f$, for each \f$ i \f$=1 .. \c m_numRefValues.
  std::vector<double> m_rateValues;
};

//*****************************************************
// Quadratic 1D->1D class
//*****************************************************
/*! \class Quadratic1D1DFunction
    \brief Class for one-dimensional quadratic functions.

    This class implements quadratic one-dimensional functions.
    A quadratic function, in mathematics, is a polynomial function of the form
    \f[ f(x)=ax^2+bx+c,\quad a \ne 0.\f]
    The graph of a quadratic function is a parabola whose axis of symmetry is parallel
    to the y-axis.
*/
class Quadratic1D1DFunction : public Base1D1DFunction {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Constructs a quadratic function in the interval <c>[minDomainValue,maxDomainValue]</c>,
   * of the type: \f$ f(x)=ax^2+bx+c \f$.*/
  Quadratic1D1DFunction(double minDomainValue,
                               double maxDomainValue,
                               double a,
                               double b,
                               double c);
  //! Destructor.
  ~Quadratic1D1DFunction();

    //! @name Mathematical  methods
  //@{
  //! Returns the value of the quadratic function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it evaluates the function at such point, namely:
   * <c>imageValue = a*domainValue^2 + b*domainValue + c </c>.*/
  double value(double domainValue) const;

  //! Returns the value of the derivative of the quadratic function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it returns the value of the derivative at such point (namely
   * <c>2*a*domainValue + b </c>. */
  double deriv(double domainValue) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  //! 'Quadratic' coefficient of the quadratic function; \f$ a \f$ in \f$ f(x)=ax^2+bx+c \f$.
  double m_a;

  //! 'Linear' coefficient of the quadratic function; \f$ b \f$ in \f$ f(x)=ax^2+bx+c \f$.
  double m_b;

  //! 'Free' coefficient of the quadratic function; \f$ c \f$ in \f$ f(x)=ax^2+bx+c \f$.
  double m_c;
};

//*****************************************************
// Sampled 1D->1D class
//*****************************************************
/*! \class Sampled1D1DFunction
    \brief Class for one-dimensional sampled functions.

    This class implements sampled one-dimensional functions.
    A sample function is one whose values are known only at discrete values of the independent
    variable; i.e., its values are only known at grid points. There are several QUESO classes
    which handle grids, all derived from BaseOneDGrid.
    Sampled functions are usually defined by arrays. */

class Sampled1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor. It should not be called by the user.
  Sampled1D1DFunction();

  //! Constructor.
  /*! When calling this constructor, the user provides the values of the independent variable
   *(\c domainValues) and their respective values in the image set (independent variable,
   \c imageValue)*/
  Sampled1D1DFunction(const std::vector<double>& domainValues,
                             const std::vector<double>& imageValues);

  //! Destructor.
  virtual ~Sampled1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the sampled function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case,  it looks for the image value associated to the domain value
   * passed to the function. If there isn't any, it calculates a linear approximation for the
   * image value of \c domainValue, considering its neighbors points in the domain.*/
  virtual double       value(double domainValue) const;

  //! <b>Bogus</b>: Derivative of the function.
  /*! Derivatives are not defined over sampled functions! Thus, this function simply checks if
   * point \c domainValue belongs to the domain of \c this function, and in affirmative case,
   * it returns 0.*/
  double               deriv(double domainValue) const;

  //! Array of the domain values (values of the independent variable).
  const   std::vector<double>& domainValues() const;

  //! Array of the image values (values of the dependent variable).
  const   std::vector<double>& imageValues () const;

  //! Checks whether the domain value \c domainValue matches exactly one of the values in the function domain \c domainValues().
  bool    domainValueMatchesExactly(double domainValue) const;

  //! Sets the values of the independent (\c domainValues) and dependent (\c imageValues) variables of this sampled function.
  void    set(const std::vector<double>& domainValues,
	      const std::vector<double>& imageValues);
  //@}

  //! @name I/O methods
  //@{
  //! Prints the values of the function in Matlab/Octave format.
  virtual void                 printForMatlab(const BaseEnvironment& env, std::ofstream& ofsvar, const std::string& prefixName) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  //! Array of the values in the domain of the function (values of the independent variable).
  std::vector<double> m_domainValues;

  //! Array of the values in the image of the function (values of the dependent variable).
  std::vector<double> m_imageValues;
};

//*****************************************************
// 'ScalarTimesFunc' 1D->1D class
//*****************************************************
/*! \class ScalarTimesFunc1D1DFunction
    \brief Class for multiplication of a one-dimensional function by a scalar.

    This class implements the multiplication of a generic one-dimensional function
    (implemented through any class derived from Base1D1DFunction) by a given
    scalar.
*/

class ScalarTimesFunc1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  ScalarTimesFunc1D1DFunction(double scalar,
                                     const Base1D1DFunction& func);

  //! Destructor.
  ~ScalarTimesFunc1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the (one-dimensional) function at point \c domainValue, already multiplied by the scalar (given during construction).
  /*! This function calls the value() method of the given function and multiplies its return
   * value by the scalar passed to the constructor. */
  double value(double domainValue) const;

  //! TODO: Returns the value of the derivative of the function multiplied by the given scalar at point \c domainValue.
  /*! \todo Please, implement me! */
  double deriv(double domainValue) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  double m_scalar;
  const Base1D1DFunction& m_func;
};

//*****************************************************
// 'FuncTimesFunc' 1D->1D class
//*****************************************************
/*! \class FuncTimesFunc1D1DFunction
    \brief Class for multiplication of a one-dimensional function by another.

    This class implements the multiplication of a generic one-dimensional function
    by another (both functions may be of any type of functions derived from
    Base1D1DFunction).
*/
class FuncTimesFunc1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  FuncTimesFunc1D1DFunction(const Base1D1DFunction& func1,
                                   const Base1D1DFunction& func2);

  //! Destructor.
  ~FuncTimesFunc1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the multiplication of a function \c func1 by another function \c func2 at the point \c domainValue, i.e. returnValue = <c>func1(domainValue)*func2(domainValue) </c>.
  /*! This function calls the value() method of the both functions and multiplies their return
   * value with one another.*/
  double value(double domainValue) const;

  //! TODO: Returns the value of the derivative of the function \c func1 by another function \c func2 at the point \c domainValue.
  /*! \todo Please, implement me! \n
   * This method objective requires clarification: are we calculating either (func1*func2)' or func1'*func2'?   */
  double deriv(double domainValue) const;
  //@}

protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  const Base1D1DFunction& m_func1;
  const Base1D1DFunction& m_func2;
};

//*****************************************************
// 'FuncPlusFunc' 1D->1D class
//*****************************************************
/*! \class FuncPlusFunc1D1DFunction
    \brief Class for addition of a one-dimensional function with another.

    This class implements the addition of two generic one-dimensional functions
    (both functions may be of any type of functions derived from Base1D1DFunction).
*/
class FuncPlusFunc1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  FuncPlusFunc1D1DFunction(const Base1D1DFunction& func1,
                                  const Base1D1DFunction& func2);

  //! Destructor.
  ~FuncPlusFunc1D1DFunction();
//@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the addition of function \c func1 and \c func2 evaluated at the point \c domainValue, i.e. returnValue = <c>func1(domainValue)+func2(domainValue) </c>.
  /*! This function calls the value() method of the both functions and adds their return
   * value with one another.*/
  double value(double domainValue) const;

  //! TODO: Returns the value of the derivative of the addition of two functions.
  /*! \todo Please, implement me! \n*/
  double deriv(double domainValue) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  const Base1D1DFunction& m_func1;
  const Base1D1DFunction& m_func2;
};
//*****************************************************
// Lagrange Polynomial 1D->1D class
//*****************************************************
/*! \class LagrangePolynomial1D1DFunction
    \brief Class for one-dimensional Lagrange polynomials.
*/

/*! The Lagrange interpolating polynomial of a one-dimensional function \f$ f(x) = y \f$
 * is the polynomial \f$ P(x) \f$ of degree \f$ \leq n-1 \f$ that passes through the
 * \f$ n \f$ points \f$ (x_1,y_1=f(x_1)),  (x_2,y_2=f(x_2)), ..., (x_n,y_n=f(x_n))\f$,
 * and is given by:
 * \f[ P(x)=\sum_{j=1}^n \prod_{k=1; k\not=j}^n \frac{x-x_k}{x_j-x_k}.\f]
 *
 * Written explicitly,
 * \f[ P(x)=\frac{(x-x_2)(x-x_3)...(x-x_n)}{(x_1-x_2)(x_1-x_3)...(x_1-x_n)}y_1+
            \frac{(x-x_1)(x-x_3)...(x-x_n)}{(x_2-x_1)(x_2-x_3)...(x_2-x_n)}y_2+...+
            \frac{(x-x_1)(x-x_2)...(x-x_(n-1))}{(x_n-x_1)(x_n-x_2)...(x_n-x_(n-1))}y_n.\f]

 * In this class, the array <c> std::vector<double>& positionValues </c> stores the points
 * \f$ x_1, x_2, ... x_n \f$  and the array <c> std::vector<double>* functionValues </c>
 * stores the points \f$ y_1, y_2, ... y_n \f$ of the Lagrange polynomial.
 *
 * \see Archer, Branden and Weisstein, Eric W. "Lagrange Interpolating Polynomial." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html.*/

class LagrangePolynomial1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  LagrangePolynomial1D1DFunction(const std::vector<double>& positionValues,
                                        const std::vector<double>* functionValues);

  //! Destructor.
  ~LagrangePolynomial1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the Lagrange polynomial at point \c domainValue.
  double value(double domainValue) const;

  //! TODO: Returns the value of the derivative of the Lagrange polynomial at point \c domainValue.
  /*! \todo This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it returns the value of the derivative at such point. */
  double deriv(double domainValue) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  std::vector<double> m_positionValues;
  std::vector<double> m_functionValues;
};

//*****************************************************
// Lagrange Basis 1D->1D class
//*****************************************************
/*! \class LagrangeBasis1D1DFunction
 *  \brief Class for Lagrange polynomial basis.
 *
 * Given a set of \f$ k+1 \f$ data points \f$(x_0, y_0),\ldots,(x_j, y_j),\ldots,(x_k, y_k)\f$
 * where no two \f$ x_j \f$ are the same, the interpolation polynomial in the Lagrange form is
 * a linear combination
 * \f[ L(x) = \sum_{j=0}^{k} y_j \ell_j(x) \f]
 * of Lagrange basis polynomials
 * \f[ \ell_j(x) = \prod_{0\le m\le k;\, m\neq j}
       \frac{x-x_m}{x_j-x_m} = \frac{(x-x_0)}{(x_j-x_0)} \cdots
       \frac{(x-x_{j-1})}{(x_j-x_{j-1})} \frac{(x-x_{j+1})}{(x_j-x_{j+1})}
       \cdots \frac{(x-x_k)}{(x_j-x_k)},\f]
  * where \f$ 0\le j\le k \f$.\n
  *
  * This class implements the one-dimensional function (Lagrange basis) \f$ \ell_j(x) \f$.
  * In this class, the array <c> std::vector<double>& positionValues </c> stores the points
  * \f$ x_1, x_2, ... x_n \f$  and the  index \f$ j \f$ is stored in \c basisIndex. */

class LagrangeBasis1D1DFunction : public Base1D1DFunction {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  LagrangeBasis1D1DFunction(const std::vector<double>& positionValues,
                                   unsigned int basisIndex);

  //! Destructor.
  ~LagrangeBasis1D1DFunction();
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the value of the Lagrange basis at point \c domainValue.
  double value(double domainValue) const;

  //! TODO: Returns the value of the derivative of the Lagrange basis at point \c domainValue.
  /*! \todo This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it returns the value of the derivative at such point. */
  double deriv(double domainValue) const;
  //@}
protected:
  using Base1D1DFunction::m_minDomainValue;
  using Base1D1DFunction::m_maxDomainValue;

  std::vector<double> m_positionValues;
  unsigned int        m_basisIndex;
};

//----------------------------------------------------------------------
//! Calculates the integral of a 2D Gaussian KDE.
/*! This function performs numerical integration (via Hermite-Gauss quadrature,
 * see GaussianHermite1DQuadrature), of a bi-variate Gaussian KDE (refer
 * to ScalarSequence::subGaussian1dKde(), ArrayOfSequences::gaussianKDE(),
 * or SequenceOfVectors::subGaussian1dKde(). */

template <class T>
double SubF1F2Gaussian2dKdeIntegral(const ScalarSequence<T>& scalarSeq1,
    const ScalarSequence<T>& scalarSeq2, unsigned int initialPos,
    double scaleValue1, double scaleValue2, const Base1D1DFunction& func1,
    const Base1D1DFunction& func2, unsigned int quadratureOrder);

}  // End namespace QUESO

#endif // UQ_1D_1D_FUNCTION_H
