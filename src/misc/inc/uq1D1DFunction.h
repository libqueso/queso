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

#ifndef __UQ_1D_1D_FUNCTION_H__
#define __UQ_1D_1D_FUNCTION_H__

#include <uqScalarSequence.h>
#include <uq1DQuadrature.h>
#include <uqEnvironment.h>
#include <uqDefines.h>
#include <vector>
#include <math.h>
#include <fstream>

/*! \file uqBase1D1DFunction.h
    \brief One-dimension function class.
*/



//*****************************************************
// Base 1D->1D class
//*****************************************************

/*! \class uqBase1D1DFunctionClass
    \brief Class for one-dimensional functions. 
    
    Base class for one-dimensional functions. 
*/
class
uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqBase1D1DFunctionClass(double minDomainValue,
			  double maxDomainValue);
  
  //! Destructor.
  virtual ~uqBase1D1DFunctionClass();
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
  virtual  double multiplyAndIntegrate(const uqBase1D1DFunctionClass& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const;
  //@}
protected:
  double m_minDomainValue;
  double m_maxDomainValue;
};

//*****************************************************
// Generic 1D->1D class
//*****************************************************
/*! \class uqGeneric1D1DFunctionClass
    \brief Class for generic one-dimensional functions. 
    
    This class implements generic one-dimensional functions. 
*/
class uqGeneric1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqGeneric1D1DFunctionClass(double minDomainValue,
                             double maxDomainValue,
                             double (*valueRoutinePtr)(double domainValue, const void* routinesDataPtr),
                             double (*derivRoutinePtr)(double domainValue, const void* routinesDataPtr),
                             const void* routinesDataPtr);
  //! Destructor.
  ~uqGeneric1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  double (*m_valueRoutinePtr)(double domainValue, const void* routinesDataPtr);
  double (*m_derivRoutinePtr)(double domainValue, const void* routinesDataPtr);
  const void* m_routinesDataPtr;
};

//*****************************************************
// Constant 1D->1D class
//*****************************************************
/*! \class uqConstant1D1DFunctionClass
    \brief Class for constant one-dimensional functions. 
    
    This class implements constant one-dimensional functions. 
*/
class uqConstant1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqConstant1D1DFunctionClass(double minDomainValue,
                              double maxDomainValue,
                              double constantValue);
  //! Destructor.
  ~uqConstant1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  double m_constantValue;
};

//*****************************************************
// Linear 1D->1D class
//*****************************************************
/*! \class uqLinear1D1DFunctionClass
    \brief Class for linear one-dimensional functions. 
    
    This class implements linear one-dimensional functions. 
    A common linear function is the the equation of a line: \f[f(x) = y_1 + m (x - x_1), \f] 
    which is a linear function with slope (or rate) \f$ m \f$ and passing through the point 
    \f$(x_1,y_1)\f$.
    
    This class implements a linear function such as the one describe above where the rate is
    given by \c rateValue, and the point is <c>(referenceDomainValue,referenceImageValue)</c>.
*/
class uqLinear1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Constructs a linear function of constant rate of change \c rateValue which passes through
   *  the point (x,y)=<c>(referenceDomainValue,referenceImageValue)</c>.*/
  uqLinear1D1DFunctionClass(double minDomainValue,
                            double maxDomainValue,
                            double referenceDomainValue,
                            double referenceImageValue,
                            double rateValue);
  //! Destructor.
  ~uqLinear1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

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
/*! \class uqPiecewiseLinear1D1DFunctionClass
    \brief Class for piecewise-linear one-dimensional functions. 
    
    This class implements piecewise-linear one-dimensional functions. 
*/
class uqPiecewiseLinear1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqPiecewiseLinear1D1DFunctionClass(double                     minDomainValue,
                                     double                     maxDomainValue,
                                     const std::vector<double>& referenceDomainValues,
                                     double                     referenceImageValue0,
                                     const std::vector<double>& rateValues);
  //! Destructor.
  ~uqPiecewiseLinear1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

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
/*! \class uqQuadratic1D1DFunctionClass
    \brief Class for generic one-dimensional functions. 
    
    This class implements quadratic one-dimensional functions. 
    A quadratic function, in mathematics, is a polynomial function of the form
    \f[ f(x)=ax^2+bx+c,\quad a \ne 0.\f]
    The graph of a quadratic function is a parabola whose axis of symmetry is parallel 
    to the y-axis.
*/
class uqQuadratic1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
    //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Constructs a quadratic function in the interval <c>[minDomainValue,maxDomainValue]</c>,
   * of the type: \f$ f(x)=ax^2+bx+c \f$.*/
  uqQuadratic1D1DFunctionClass(double minDomainValue,
                               double maxDomainValue,
                               double a,
                               double b,
                               double c);
  //! Destructor.
  ~uqQuadratic1D1DFunctionClass();

    //! @name Mathematical  methods
  //@{
  //! Returns the value of the quadratic function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it evaluates the function at such point, namely: 
   * <c>imageValue = a*domainValue^2 + b*domainValue + c </c>\f$.*/
  double value(double domainValue) const;
  
  //! Returns the value of the derivative of the quadratic function at point \c domainValue.
  /*! This function checks if point \c domainValue belongs to the domain of \c this function,
   * and in affirmative case, it returns the value of the derivative at such point (namely 
   * <c>2*a*domainValue + b </c>. */
  double deriv(double domainValue) const;
  //@}
protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  //! 'Quadratic' coefficient of the quadractic function; \f$ a \f$ in \f$ f(x)=ax^2+bx+c \f$.
  double m_a;
  
  //! 'Linear' coefficient of the quadractic function; \f$ b \f$ in \f$ f(x)=ax^2+bx+c \f$.
  double m_b;
  
  //! 'Free' coefficient of the quadractic function; \f$ c \f$ in \f$ f(x)=ax^2+bx+c \f$.
  double m_c;
};

//*****************************************************
// Sampled 1D->1D class
//*****************************************************
/*! \class uqSampled1D1DFunctionClass
    \brief Class for one-dimensional sampled functions. 
    
    This class implements sampled one-dimensional functions. 
    A sample function is one whose values are known only at discrete values of the independent
    variable; i.e., its values are only known at grid points. There are several QUESO classes 
    which handle grids, all derived from uqBaseOneDGridClass.  
    Sampled functions are usually defined by arrays. */

class uqSampled1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor. It should not be called by the user.
  uqSampled1D1DFunctionClass();
  
  //! Constructor. 
  /*! When calling this constructor, the user provides the values of the independent variable 
   *(\c domainValues) and their respective values in the image set (independent variable, 
   \c imageValue)*/
  uqSampled1D1DFunctionClass(const std::vector<double>& domainValues,
                             const std::vector<double>& imageValues);
  
  //! Destructor.
  virtual ~uqSampled1D1DFunctionClass();
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
  virtual void                 printForMatlab(const uqBaseEnvironmentClass& env, std::ofstream& ofsvar, const std::string& prefixName) const;
  //@}
protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  //! Array of the values in the domain of the function (values of the independent variable).
  std::vector<double> m_domainValues;
  
  //! Array of the values in the image of the function (values of the dependent variable).
  std::vector<double> m_imageValues;
};

//*****************************************************
// Delta Set 1D->1D class
//*****************************************************
/*!\class uqDeltaSet1D1DFunctionClass
 * \brief A class that implements the Delta function.
 *
 * In mathematics, the Dirac delta function, or \f$ \delta \f$ function, is 
 * (informally) a generalized function on the real number line that is zero 
 * everywhere except at zero, with an integral of one over the entire real line.
 * The delta function has the fundamental property that:
 *  \f[ \int_{-\infty}^{\infty} f(x) \delta(x-a)\, dx=f(a) \f]
 * and, in fact,
 *  \f[ \int_{a-\epsilon}^{a+\epsilon} f(x) \delta(x-a)\,dx=f(a)\f] 
 * for  \f$ \epsilon > 0 \f$.
 */
class uqDeltaSet1D1DFunctionClass : public uqBase1D1DFunctionClass { // july2011
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqDeltaSet1D1DFunctionClass(const std::vector<double>& domainValues,
                              double                     domainMin,
                              double                     domainMax,
                              const std::vector<double>& imageValues);
  //! Destructor.
  ~uqDeltaSet1D1DFunctionClass();
  //@}
  
  //! @name Mathematical  methods
  //@{
  //! Returns the value of the Delta function point \c domainValue.
  double value               (double domainValue) const;
  
  //! TODO: Returns the value of the derivative of the Delta function  (which is always zero) point \c domainValue.
  double deriv               (double domainValue) const;
  
  //! TODO: Multiplies the Delta function with another function \c func and integrates it.
  /*! \todo This method calculates numerically the integral: 
   * \f$ \int f(x) \delta(x)\, dx */
  double multiplyAndIntegrate(const uqBase1D1DFunctionClass& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const; // july2011
  //@}
protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  std::vector<double> m_domainValues;
  std::vector<double> m_imageValues;

  mutable int                 m_lastInputK;
  mutable unsigned int        m_lastQuadratureOrder;
  mutable std::vector<double> m_level_0;
  mutable std::vector<double> m_level_m1;
};

//*****************************************************
// Gaussian 1d Kde 1D->1D class
//*****************************************************
/*! \class uqGaussian1dKde1D1DFunctionClass
    \brief Class for one-dimensional kernel-density estimator (functions) using Gaussian kernels. 
    
    This class handles the representation of a kernel-density estimator function that 
    uses Gaussian kernels.
*/

class uqGaussian1dKde1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! \p chain: contains the data points  from which the density estimate is constructed,
   * \p chainMin and \p chainMax - define the interval <c>[chainMin,chainMax]</c> on which 
   * the density estimate is constructed.*/
  uqGaussian1dKde1D1DFunctionClass(const uqScalarSequenceClass<double>* chain,
                                   double chainMin,
                                   double chainMax,
                                   double gaussian1DScale);
  //! Destructor.
  ~uqGaussian1dKde1D1DFunctionClass();
  //@}
  
  //! @name Mathematical  methods
  //@{ 
  //! Returns the value of the Gaussian KDE estimator at point \c domainValue.
  /*! Let \f$(x_1, x_2, …, x_n)\f$ be an  independent and identically distributed (iid) 
     sample drawn from some distribution with an unknown density \f$ f \f$.  Its kernel 
     density estimator is:  \f[\hat{f}_h(x) = \frac{1}{n}\sum_{i=1}^n K_h (x - x_i) = 
     \frac{1}{nh} \sum_{i=1}^n K\Big(\frac{x-x_i}{h}\Big),\f] 
     where \f$ K() \f$ is the kernel, and \f$ h > 0 \f$ is a smoothing parameter (scale)
     called the bandwidth. \n In the present function, the RHS of the equation
     above is implemented using function uqMiscGaussianDensity() as the Gaussian kernel, 
     and the scale parameter is stored in \c m_gaussian1DScale. */
  double value(double domainValue) const;
  
  //! TODO: Returns the value of the derivative of the Gaussian KDE function.
  /*! \todo Please, implement me!*/
  double deriv(double domainValue) const;
  
  //! Multiplies \c this function with \c function, and integrates it numerically, using a scheme of order \c quadratureOrder. 
  double multiplyAndIntegrate(const uqBase1D1DFunctionClass& func, unsigned int quadratureOrder, double* resultWithMultiplicationByTAsWell) const;
  //@}
  
protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  //! Datapoints from which the density estimate is constructed,
  const uqScalarSequenceClass<double>* m_chain;
  
  //! Minimum value of the interval <c>[chainMin,chainMax]</c> on which the density estimate is constructed.
  double m_chainMin;
  
  //! Maximum value of the interval <c>[chainMin,chainMax]</c> on which the density estimate is constructed.
  double m_chainMax;
  
  //! Scale/bandwidth of the Gaussian (smoothing parameter). 
  double m_gaussian1DScale;

  mutable int                 m_lastInputK;
  mutable unsigned int        m_lastQuadratureOrder;
  mutable std::vector<double> m_level_0;
  mutable std::vector<double> m_level_m1;
};

//*****************************************************
// 'ScalarTimesFunc' 1D->1D class
//*****************************************************
/*! \class uqScalarTimesFunc1D1DFunctionClass
    \brief Class for multiplication of a one-dimensional function by a scalar. 
    
    This class implements the multiplication of a generic one-dimensional function 
    (implemented through any class derived from uqBase1D1DFunctionClass) by a given
    scalar. 
*/

class uqScalarTimesFunc1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqScalarTimesFunc1D1DFunctionClass(double scalar,
                                     const uqBase1D1DFunctionClass& func);
 
  //! Destructor.
  ~uqScalarTimesFunc1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  double m_scalar;
  const uqBase1D1DFunctionClass& m_func;
};

//*****************************************************
// 'FuncTimesFunc' 1D->1D class
//*****************************************************
/*! \class uqFuncTimesFunc1D1DFunctionClass
    \brief Class for multiplication of a one-dimensional function by another. 
    
    This class implements the multiplication of a generic one-dimensional function 
    by another (both functions may be of any type of functions derived from 
    uqBase1D1DFunctionClass). 
*/
class uqFuncTimesFunc1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqFuncTimesFunc1D1DFunctionClass(const uqBase1D1DFunctionClass& func1,
                                   const uqBase1D1DFunctionClass& func2);
 
  //! Destructor.
  ~uqFuncTimesFunc1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqBase1D1DFunctionClass& m_func1;
  const uqBase1D1DFunctionClass& m_func2;
};

//*****************************************************
// 'FuncPlusFunc' 1D->1D class
//*****************************************************
/*! \class uqFuncPlusFunc1D1DFunctionClass
    \brief Class for addition of a one-dimensional function with another. 
    
    This class implements the addition of two generic one-dimensional functions 
    (both functions may be of any type of functions derived from uqBase1D1DFunctionClass). 
*/
class uqFuncPlusFunc1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqFuncPlusFunc1D1DFunctionClass(const uqBase1D1DFunctionClass& func1,
                                  const uqBase1D1DFunctionClass& func2);
  
  //! Destructor.
  ~uqFuncPlusFunc1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqBase1D1DFunctionClass& m_func1;
  const uqBase1D1DFunctionClass& m_func2;
};

//*****************************************************
// 'QuadDenominator' 1D->1D class
//*****************************************************
/*! \class uqQuadDenominator1D1DFunctionClass
    \brief Class for multiplication of a one-dimensional function with itself. 
    
    Given a generic one-dimensional function \f$ f=f(x)\f$ (implemented through 
    any class derived from uqBase1D1DFunctionClass), this class implements another
    one-dimensional function \f$ g=g(x) \f$ of the form:
    \f[ g(x) = f^2(x). \f]
    In other words, this class implements the squared of a generic one-dimensional 
    function.*/

class uqQuadDenominator1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqQuadDenominator1D1DFunctionClass(const uqBase1D1DFunctionClass& func);
  
  //! Destructor.
  ~uqQuadDenominator1D1DFunctionClass();
  //@}
  
  //! @name Mathematical  methods
  //@{
  //! Returns square of the value of \c this function evaluated at point \c domainValue.
  /*! This function calls the value() method of the given function and returns its squared value. */
  double value(double domainValue) const;
  
  //! TODO: Returns the value of the derivative \f$ g' \f$ of the function, \f$ g(x) = f^2(x)\f$.
  /*! \todo Please, implement me! \n  
   * This method objective requires clarification: 
   * are we calculating either (f^2)' or (f')^2 */  
  double deriv(double domainValue) const;
  
  //! Access to the constitutive function passed to the constructor; internally stored in \c m_func.
  const uqBase1D1DFunctionClass& constitutiveFunction() const;
  //@}
protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqBase1D1DFunctionClass& m_func;
};

//*****************************************************
// 'QuadNumerator' 1D->1D class
//*****************************************************
/*! \class uqQuadNumerator1D1DFunctionClass
    \brief Class for multiplication of a one-dimensional function with both itself and the independent varible. 
    
    Given a generic one-dimensional function \f$ f=f(x)\f$, this class implements another
    one-dimensional function \f$ g=g(x) \f$ of the form:
    \f[ g(x) = x f^2(x). \f]
  */
class uqQuadNumerator1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqQuadNumerator1D1DFunctionClass(const uqBase1D1DFunctionClass& func);

  //! Destructor.
  ~uqQuadNumerator1D1DFunctionClass();
  //@}
  
  //! @name Mathematical  methods
  //@{
  //! Returns the squared value of the (one-dimensional) squared function at point \c domainValue multiplied by \c domainValue.
  /*! This function calls the value() method of the given function and returns its squared value multiplied by \c domainValue. */
  double value(double domainValue) const;
  
  //! TODO: Returns the value of the derivative of the function \f$ g(x) = x f^2(x)\f$.
  /*! \todo Please, implement me! 
   * This method objective requires clarification: 
   * are we calculating either (xf^2)' or x(f')^2 */    
  double deriv(double domainValue) const;
  //@}
protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  const uqBase1D1DFunctionClass& m_func;
};

//*****************************************************
// 'QuadCRecursion' 1D->1D class
//*****************************************************
class uqQuadCRecursion1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  uqQuadCRecursion1D1DFunctionClass(int k);
  uqQuadCRecursion1D1DFunctionClass(int                        k,
                                    const std::vector<double>& alpha,
                                    const std::vector<double>& beta);
 ~uqQuadCRecursion1D1DFunctionClass();

        double               value(double domainValue) const;
        double               deriv(double domainValue) const;
        int                  k    ()                   const;
  const std::vector<double>& alpha()                   const;
  const std::vector<double>& beta ()                   const;

protected:
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;
  int                 m_k;
  std::vector<double> m_alpha;
  std::vector<double> m_beta;
};

//*****************************************************
// Templated isolated function
//*****************************************************
template<class V>
void
computeAlphasAndBetasWithRiemann(
  const uqBase1D1DFunctionClass& function1D1D,
  unsigned int numIntegrationPoints,
  V&           alpha,
  V&           beta)
{
  unsigned int n = alpha.sizeLocal();
  UQ_FATAL_TEST_MACRO((n < 1),
                      UQ_UNAVAILABLE_RANK,
                      "computeAlphasAndBetasWithRiemann()",
                      "invalid n");

  std::vector<double> samplePositions(numIntegrationPoints,0.);
  std::vector<double> sampleValues   (numIntegrationPoints,0.);

  double delta = (function1D1D.maxDomainValue()-function1D1D.minDomainValue())/((double) numIntegrationPoints);
  for (unsigned int i = 0; i < numIntegrationPoints; ++i) {
    samplePositions[i] = function1D1D.minDomainValue() + (.5 + ((double) i))*delta;
    sampleValues   [i] = function1D1D.value(samplePositions[i]);
    UQ_FATAL_TEST_MACRO((sampleValues[i] < 0.),
                        UQ_UNAVAILABLE_RANK,
                        "computeAlphasAndBetasWithRiemann<V>()",
                        "sampleValue is negative");
  }

  std::vector<double> pi_m1(numIntegrationPoints,0.);
  std::vector<double> pi_0 (numIntegrationPoints,1.); // Yes, '1'
  std::vector<double> pi_p1(numIntegrationPoints,0.);
  double pi_pi_m1 = 0.;

  for (unsigned int k = 0; k < n; ++k) {
    double pi_pi_0   = 0.;
    double t_pi_pi_0 = 0.;
    for (unsigned int i = 0; i < numIntegrationPoints; ++i) {
      pi_pi_0   += delta                      * pi_0[i] * pi_0[i] * sampleValues[i];
      t_pi_pi_0 += delta * samplePositions[i] * pi_0[i] * pi_0[i] * sampleValues[i];
    }

    // Alpha and beta
    alpha[k] = t_pi_pi_0/pi_pi_0;
    if (k == 0) {
      beta[k] = 0.;
      for (unsigned int i = 0; i < numIntegrationPoints; ++i) {
        beta[k] += delta * sampleValues[i];
      }
    }
    else {
      beta[k] = pi_pi_0/pi_pi_m1;
    }
    UQ_FATAL_TEST_MACRO((beta[k] < 0.),
                        UQ_UNAVAILABLE_RANK,
                        "computeAlphasAndBetasWithRiemann<V>()",
                        "beta is negative");

    // Prepare for next k
    if (k < (n-1)) {
      for (unsigned int i = 0; i < numIntegrationPoints; ++i) {
        pi_p1[i] = (samplePositions[i] - alpha[k])*pi_0[i] - beta[k]*pi_m1[i];
      }
      pi_m1    = pi_0;
      pi_0     = pi_p1;
      pi_pi_m1 = pi_pi_0;
    }
  }

  return;
}

//*****************************************************
// Isolated function
//*****************************************************
void
alphaBetaCLoop( // july2011
  int                                      k,
  double                                   pi_pi_m1,
  const uqQuadCRecursion1D1DFunctionClass& pi_m1,
  const uqQuadCRecursion1D1DFunctionClass& pi_0,
  const uqBase1D1DFunctionClass&           rho,
  unsigned int                             integrationOrder,
  std::vector<double>&                     alpha,
  std::vector<double>&                     beta);

//*****************************************************
// Isolated templated function
//*****************************************************
template<class V>
void
computeAlphasAndBetasWithCQuadrature( // july2011
  const uqBase1D1DFunctionClass& function1D1D,
  unsigned int                   integrationOrder,
  V&                             alpha,
  V&                             beta)
{
  uqQuadCRecursion1D1DFunctionClass pi_m1(-1);
  uqQuadCRecursion1D1DFunctionClass pi_0 (0);

  std::vector<double> stdAlpha(alpha.sizeLocal(),0.);
  std::vector<double> stdBeta (beta.sizeLocal(), 0.);
  alphaBetaCLoop(0, // Yes, '0'
                 1.,
                 pi_m1,
                 pi_0,
                 function1D1D,
                 integrationOrder,
                 stdAlpha,
                 stdBeta);

  for (unsigned int i = 0; i < alpha.sizeLocal(); ++i) {
    alpha[i] = stdAlpha[i];
    beta [i] = stdBeta [i];
  }

  return;
}

//#include <uq1D1DFunction2.h>

//*****************************************************
// Isolated templated function
//*****************************************************
template<class V, class M>
void
computeQuadPtsAndWeights(const uqBase1D1DFunctionClass& function1D1D,
                         unsigned int                   integrationMethod,
                         unsigned int                   integrationOrder,
                         double                         exactValueToForceOnSumOfWeights,
                         V&                             quadPositions,
                         V&                             quadWeights)
{
  const uqBaseEnvironmentClass& env = quadPositions.env();
  unsigned int n = quadPositions.sizeLocal();

  UQ_FATAL_TEST_MACRO((n == 0),
                      UQ_UNAVAILABLE_RANK,
                      "computeQuadPtsAndWeights<V,M>()",
                      "invalid input vector size");

  UQ_FATAL_TEST_MACRO((quadPositions.sizeLocal() != quadWeights.sizeLocal()),
                      UQ_UNAVAILABLE_RANK,
                      "computeQuadPtsAndWeights<V,M>()",
                      "different input vector sizes");

  // Compute alphas and betas
  V alpha(quadPositions);
  V beta (quadPositions);
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>()"
                          << ": computing alphas and betas ..."
                          << std::endl;
  }

  if (integrationMethod == 0) { // Riemann
    computeAlphasAndBetasWithRiemann<V>(function1D1D,
                                        integrationOrder,
                                        alpha,
                                        beta);
  }
  else if (integrationMethod == 1) { // Quadrature
    //computeAlphasAndBetasWithQuadrature<V>   (function1D1D, // IMPORTANT
    //computeAlphasAndBetasWithPolQuadrature<V>(function1D1D, // IMPORTANT
    computeAlphasAndBetasWithCQuadrature<V>(function1D1D, // IMPORTANT
                                            integrationOrder,
                                            alpha,
                                            beta);
  }
  else {
    UQ_FATAL_TEST_MACRO(true,
                        UQ_UNAVAILABLE_RANK,
                        "computeQuadPtsAndWeights<V,M>()",
                        "invalid 'quadIntegrationMethod'");
  }

  if ((env.displayVerbosity() >= 2) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>()";
    for (unsigned int k = 0; k < n; ++k) {
      *env.subDisplayFile() << "\n alpha[" << k << "] = " << alpha[k]
                            << ", beta["   << k << "] = " << beta[k];
    }
    *env.subDisplayFile() << std::endl;
  }

  // Form mJ
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>()"
                          << ": forming J..."
                          << std::endl;
  }
  V zeroVector(quadPositions);
  zeroVector.cwSet(0.);
  M mJ(zeroVector);
  for (unsigned int k = 0; k < n; ++k) {
    mJ(k,k) = alpha[k];
    if (mJ.numRowsGlobal() > 1) {
      if (k < (mJ.numRowsGlobal()-1)) {      
        mJ(k,k+1) = sqrt(beta[k+1]);
        mJ(k+1,k) = sqrt(beta[k+1]);
      }
    }
  } // end for 'k'
  if ((env.displayVerbosity() >= 2) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>():"
                          << "\n  mJ = " << mJ
                          << std::endl;
  }

  // Compute eigen values and vectors of mJ
  if ((env.displayVerbosity() >= 3) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "In computeQuadPtsAndWeights<V,M>()"
                          << ": computing eigen stuff of J..."
                          << std::endl;
  }
  V eigenValues(quadPositions);
  M eigenVectors(zeroVector);
  mJ.eigen(eigenValues,&eigenVectors);

  // Prepare output information
  for (unsigned int k = 0; k < n; ++k) {
    quadPositions[k] = eigenValues[k];
    quadWeights  [k] = beta[0] * eigenVectors(0,k) * eigenVectors(0,k);
  } 

  double wSum = 0.;
  for (unsigned int k = 0; k < n; ++k) {
    wSum += quadWeights[k];
  } 

  if ((env.displayVerbosity() >= 2) && (env.subDisplayFile())) {
    *env.subDisplayFile() << "IMPORTANT: In computeQuadPtsAndWeights<V,M>()"
                          << ", before consulting 'exactValueToForceOnSumOfWeights'"
                          << ": exactValueToForceOnSumOfWeights = " << exactValueToForceOnSumOfWeights
                          << ", beta[0] = "                         << beta[0]
                          << ", wSum = "                            << wSum
                          << "\n  eigenValues = "                   << eigenValues
                          << "\n  eigenVectors = "                  << eigenVectors
                          << std::endl;
  }

  if (exactValueToForceOnSumOfWeights > 0.) {
    for (unsigned int k = 0; k < n; ++k) {
      quadWeights[k] *= (exactValueToForceOnSumOfWeights/wSum);
    }
  }

  return;
}

//*****************************************************
// Lagrange Polynomial 1D->1D class
//*****************************************************
/*! \class uqLagrangePolynomial1D1DFunctionClass
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

class uqLagrangePolynomial1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqLagrangePolynomial1D1DFunctionClass(const std::vector<double>& positionValues,
                                        const std::vector<double>* functionValues);
 
  //! Destructor.
  ~uqLagrangePolynomial1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  std::vector<double> m_positionValues;
  std::vector<double> m_functionValues;
};

//*****************************************************
// Lagrange Basis 1D->1D class
//*****************************************************
/*! \class uqLagrangeBasis1D1DFunctionClass
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

class uqLagrangeBasis1D1DFunctionClass : public uqBase1D1DFunctionClass {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqLagrangeBasis1D1DFunctionClass(const std::vector<double>& positionValues,
                                   unsigned int basisIndex);
 
  //! Destructor.
  ~uqLagrangeBasis1D1DFunctionClass();
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
  using uqBase1D1DFunctionClass::m_minDomainValue;
  using uqBase1D1DFunctionClass::m_maxDomainValue;

  std::vector<double> m_positionValues;
  unsigned int        m_basisIndex;
};

//----------------------------------------------------------------------
//! Calculates the integral of a 2D Gaussian KDE.
/*! This function performs numerical integration (via Hermite-Gauss quadrature, 
 * see uqGaussianHermite1DQuadratureClass), of a bi-variate Gaussian KDE (refer
 * to uqScalarSequenceClass::subGaussian1dKde(), uqArrayOfSequencesClass::gaussianKDE(), 
 * or uqSequenceOfVectorsClass::subGaussian1dKde(). */

template <class T>
double
uqSubF1F2Gaussian2dKdeIntegral(const uqScalarSequenceClass<T>& scalarSeq1,
                               const uqScalarSequenceClass<T>& scalarSeq2,
                               unsigned int                    initialPos,
                               double                          scaleValue1,
                               double                          scaleValue2,
                               const uqBase1D1DFunctionClass&  func1,
                               const uqBase1D1DFunctionClass&  func2,
                               unsigned int                    quadratureOrder)
{
  double resultValue = 0.;

  UQ_FATAL_TEST_MACRO(initialPos != 0,
                      scalarSeq1.env().worldRank(),
                      "uqSubF1F2Gaussian2dKdeIntegral()",
                      "not implemented yet for initialPos != 0");
  UQ_FATAL_TEST_MACRO(scalarSeq1.subSequenceSize() != scalarSeq2.subSequenceSize(),
                      scalarSeq1.env().worldRank(),
                      "uqSubF1F2Gaussian2dKdeIntegral()",
                      "different sizes");

  uqGaussianHermite1DQuadratureClass quadObj(0.,1.,quadratureOrder);
  const std::vector<double>& quadPositions = quadObj.positions();
  const std::vector<double>& quadWeights   = quadObj.weights  ();
  UQ_FATAL_TEST_MACRO(quadPositions.size() != quadWeights.size(),
                      UQ_UNAVAILABLE_RANK,
                      "uqSubF1F2Gaussian2dKdeIntegral()",
                      "quadObj has invalid state");

  unsigned int numQuadraturePositions = quadPositions.size();
  unsigned int dataSize = scalarSeq1.subSequenceSize();
  for (unsigned int k = 0; k < dataSize; ++k) {
    double value1 = 0.;
    double value2 = 0.;
    double x1k = scalarSeq1[k];
    double x2k = scalarSeq2[k];
    for (unsigned int j = 0; j < numQuadraturePositions; ++j) {
      value1 += func1.value(scaleValue1*quadPositions[j]+x1k)*quadWeights[j];
      value2 += func2.value(scaleValue2*quadPositions[j]+x2k)*quadWeights[j];
    }
    resultValue += value1*value2;
  }
  resultValue *= 1./(2.*M_PI)/((double) dataSize);

  return resultValue;
}
#endif // __UQ_1D_1D_FUNCTION_H__

