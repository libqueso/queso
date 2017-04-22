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

#ifndef UQ_1D_1D_QUADRATURE_H
#define UQ_1D_1D_QUADRATURE_H

#include <queso/BaseQuadrature.h>
#include <queso/Environment.h>
#include <queso/Defines.h>
#include <vector>
#include <math.h>
#include <fstream>

namespace QUESO {

/*!\file 1DQuadrature.h
 * \brief One-dimensional quadrature rules (numerical integration) class. */

//*****************************************************
// Base 1D quadrature class
//*****************************************************
/*! \class Base1DQuadrature
    \brief Base class for one-dimensional quadrature rules (numerical integration of functions).

    Base class for numerical integration via quadrature rules of one-dimensional functions.
*/
class Base1DQuadrature : public BaseQuadrature
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  Base1DQuadrature(double minDomainValue,
			  double maxDomainValue,
			  unsigned int order);

  //! Pure virtual destructor, forcing this to be an abstract object.
  virtual ~Base1DQuadrature() =0;
  //@}

  //! @name Mathematical  methods
  //@{
  //! Returns the minimum value of the domain of the (one-dimensional) function.
  double                     minDomainValue() const;

  //! Returns the maximum value of the domain of the (one-dimensional) function.
  double                     maxDomainValue() const;

  //! Returns the order of the quadrature rule.
  unsigned int               order         () const;

  //! Array of the positions for the numerical integration.
  const std::vector<double> & positions() const
  { queso_assert(!m_positions.empty());
    return m_positions; }
  //@}

protected:
  double              m_minDomainValue;
  double              m_maxDomainValue;
  unsigned int        m_order;
  std::vector<double> m_positions;
};

//*****************************************************
// Generic 1D quadrature class
//*****************************************************
/*! \class Generic1DQuadrature
    \brief Class for one-dimensional generic quadrature rules (numerical integration of functions).

    Class for generic quadrature rules for numerical integration of one-dimensional functions.
*/
class Generic1DQuadrature : public Base1DQuadrature {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  Generic1DQuadrature(double minDomainValue,
                             double maxDomainValue,
                             const std::vector<double>& positions,
                             const std::vector<double>& weights);

  //! Destructor.
  ~Generic1DQuadrature();
 //@}

};

//*****************************************************
// Uniform/Legendre 1D quadrature class
//*****************************************************
/*! \class UniformLegendre1DQuadrature
 *    \brief Class for Legendre-Gauss quadrature rule for one-dimensional functions.
 *
 * In a general Gaussian quadrature rule, an definite integral of \f$ f(x)\f$ is first
 * approximated over the interval [-1,1] by a polynomial approximable function \f$ g(x)\f$
 * and a known weighting function \f$ W(x)\f$:
 * \f[\int_{-1}^1 f(x) \, dx = \int_{-1}^1 W(x) g(x) \, dx\f]
 * Those are then approximated by a sum of function values at specified points \f$ x_i \f$
 * multiplied by some weights \f$ w_i \f$:
 * \f[ \int_{-1}^1 W(x) g(x) \, dx \approx \sum_{i=1}^n w_i g(x_i) \f]
 * In the case of Gauss-Legendre quadrature, the weighting function \f$ W(x) = 1 \f$,
 * so we can approximate an integral of \f$ f(x) \f$ with:
 * \f[ \int_{-1}^1 f(x)\,dx \approx \sum_{i=1}^n w_i f(x_i) \f]
 * The abscissas for quadrature order \f$ n \f$ are given by the roots of the Legendre
 * polynomials \f$ P_n(x)\f$, which occur symmetrically about 0. The weights are
 *   \f[ w_i = \frac{2}{(1-x_i^2)[P'_n(x_i)]^2}=\frac{2(1-x_i^2)}{(n+1)^2[P_{n+1}(x_i)]^2} \f]
 *
 * Several authors give a table of abscissas and weights:
 *
<table border="2">
<tr>
<th>\f$ n \f$</th><th>\f$ x_i \f$</th><th>\f$ w_i \f$</th>
</tr>
<tr><td> 2 </td><td> \f$ \pm \frac{1}{3}\sqrt{3} \f$ 		      </td><td> \f$ 1 \f$ </td></tr>
<tr><td> 3 </td><td> \f$ 0	 \f$ 				      </td><td> \f$ \frac{8}{9} \f$ </td></tr>
<tr><td>   </td><td> \f$ \pm \frac{1}{5} \sqrt{15} \f$ 	      </td><td> \f$ \frac{5}{9} \f$ </td></tr>
<tr><td> 4 </td><td> \f$ \pm \frac{1}{35}\sqrt{525-70\sqrt{30}} \f$ </td><td> \f$ \frac{1}{36}(18+\sqrt{30})\f$ </td></tr>
<tr><td>   </td><td> \f$ \pm \frac{1}{35}\sqrt{525+70\sqrt{30}} \f$ </td><td> \f$ \frac{1}{36}(18-\sqrt{30})\f$ </td></tr>
<tr><td> 5 </td><td> \f$ 0 \f$ 					      </td><td> \f$ \frac{128}{225}\f$ </td></tr>
<tr><td>   </td><td> \f$ \pm \frac{1}{21}\sqrt{245-14\sqrt{70}} \f$ </td><td> \f$ \frac{1}{900}(322+13\sqrt{70})\f$ </td></tr>
<tr><td>   </td><td> \f$ \pm \frac{1}{21}\sqrt{245+14\sqrt{70}} \f$ </td><td> \f$ \frac{1}{900}(322-13\sqrt{70})\f$ </td></tr>
</table>
 *
 * \see Weisstein, Eric W. "Legendre-Gauss Quadrature." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/Legendre-GaussQuadrature.html.*/

class UniformLegendre1DQuadrature : public Base1DQuadrature {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Constructs a Gaussian-Legendre quadrature of order \c order, in the interval
   * <c>[minDomainValue,maxDomainValue]</c>. Valid values for the order of the
   * quadrature rule are: 1-7, 10-12, 16. This method scales the abscissas (positions)
   * of the quadrature from the interval [-1,1] to <c>[minDomainValue,maxDomainValue]</c>,
   * and the parameter \c densityIsNormalized determines whether the weights should be
   * scaled as well. */
  UniformLegendre1DQuadrature(double       minDomainValue,
                                     double       maxDomainValue,
                                     unsigned int order,
                                     bool         densityIsNormalized);
  //! Destructor.
  ~UniformLegendre1DQuadrature();
  //@}

};

//*****************************************************
// Gaussian/Hermite 1D quadrature class
//*****************************************************
/*! \class GaussianHermite1DQuadrature
 *  \brief Class for Hermite-Gauss quadrature rule for one-dimensional functions.
 *
 * Hermite-Gauss quadrature, also called Hermite quadrature, is a Gaussian quadrature
 * over the interval \f$(-\infty,\infty)\f$ with weighting function \f$ W(x)=e^{-x^2}\f$.\n
 * The abscissas for quadrature order \f$ n \f$ are given by the roots \f$ x_i \f$ of
 * the Hermite polynomials \f$ H_n(x)\f$, which occur symmetrically about 0.\n
 *
 * The abscissas and weights can be computed analytically for small \f$ n \f$:
 *
 * <table border="2">
<tr>
<th>\f$ n \f$</th><th>\f$ x_i \f$</th><th>\f$ w_i \f$</th>
</tr>
<tr><td> 2 </td><td>\f$\pm \frac{1}{2}\sqrt{2} \f$ </td><td> \f$ \frac{1}{2}\sqrt{\pi} \f$ </td></tr>
<tr><td> 3 </td><td>\f$ 0 \f$			   </td><td> \f$ \frac{2}{3}\sqrt{\pi} \f$ </td></tr>
<tr><td>   </td><td>\f$\pm \frac{1}{2}\sqrt{6} \f$ </td><td> \f$ \frac{1}{6}\sqrt{\pi} \f$ </td></tr>
<tr><td> 4 </td><td>\f$\pm \sqrt{\frac{3-\sqrt{6}}{2}} \f$ </td><td> \f$ \frac{\sqrt{\pi}}{4(3-\sqrt{6})} \f$ </td></tr>
<tr><td>   </td><td>\f$\pm \sqrt{\frac{3-\sqrt{6}}{2}} \f$ </td><td> \f$ \frac{\sqrt{\pi}}{4(3+\sqrt{6})} \f$ </td></tr>
</table>
 *  \see Weisstein, Eric W. "Hermite-Gauss Quadrature." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/Hermite-GaussQuadrature.html.*/

class GaussianHermite1DQuadrature : public Base1DQuadrature {
public:
   //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Constructs a Gaussian-Hermite quadrature of order \c order.
   * Valid values for the order of the  quadrature rule are: 1-9, 19.
   * \todo: Prepare the code to include both parameters \c mean and \c stddev. */
  GaussianHermite1DQuadrature(double       mean,
                                     double       stddev,
                                     unsigned int order);
  //! Destructor.
  ~GaussianHermite1DQuadrature();
  //@}

protected:

  double m_mean;
  double m_stddev;
};

//*****************************************************
// Wigner/Chebyshev1st 1D quadrature class
//*****************************************************
/*! \class WignerInverseChebyshev1st1DQuadrature
 *  \brief Class for first type Chebyshev-Gauss quadrature rule for one-dimensional functions.
 *
 * Chebyshev-Gauss quadrature, also called Chebyshev Type 1 quadrature, is a Gaussian
 * quadrature over the interval [-1,1] with weighting function \f$ W(x)=\frac{1}{\sqrt{1-x^2}}\f$.\n
 * The abscissas for quadrature order \f$ n \f$ are given by the roots of the Chebyshev
 * polynomial of the first kind \f$ T_n(x) \f$, which occur symmetrically about 0.\n
 *
 * The abscissas are given explicitly by \f$ x_i=\cos[\frac{(2i-1)\pi}{2n}]\f$ and the
 * weights are \f$ w_i=\frac{\pi}{n}. \f$
 *
 * \see Weisstein, Eric W. "Chebyshev-Gauss Quadrature." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/Chebyshev-GaussQuadrature.html.
 * \see http://en.wikipedia.org/wiki/Chebyshev-Gauss_quadrature. */

class WignerInverseChebyshev1st1DQuadrature : public Base1DQuadrature {
public:
   //! @name Constructor/Destructor methods
  //@{
  //! TODO: Default constructor.
  /*! \todo Constructs a Gaussian-Chebyshev quadrature (of first type) of order
   * \c order, in the interval <c>[minDomainValue,maxDomainValue]</c>. This method
   * scales the the abscissas (positions) of the quadrature from the interval [-1,1]
   * to <c>[minDomainValue,maxDomainValue]</c>. */
  WignerInverseChebyshev1st1DQuadrature(double       minDomainValue,
                                               double       maxDomainValue,
                                               unsigned int order);
  //! Destructor.
  ~WignerInverseChebyshev1st1DQuadrature();
  //@}

};

//*****************************************************
// Wigner/Chebyshev2nd 1D quadrature class
//*****************************************************
/*! \class WignerChebyshev2nd1DQuadrature
 *  \brief Class for second type Chebyshev-Gauss quadrature rule for one-dimensional functions.
 *
 * Chebyshev-Gauss quadrature, also called Chebyshev Type 2 quadrature, is a Gaussian
 * quadrature over the interval [-1,1] with weighting function \f$ W(x)=\sqrt{1-x^2}\f$.\n
 * The abscissas for quadrature order \f$ n \f$ are given by the roots of the Chebyshev
 * polynomial of the \b second kind \f$ U_n(x) \f$, which occur symmetrically about 0.\n
 *
 * The abscissas are given explicitly by \f$ x_i=\cos[\frac{i\pi}{n+1}].\f$
 * and all the weights are \f$ w_i=\frac{\pi}{n+1}\sin^2[\frac{i\pi}{n+1}]. \f$
 *
 * \see  Weisstein, Eric W. "Gaussian Quadrature." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/GaussianQuadrature.html.
 * \see http://en.wikipedia.org/wiki/Chebyshev-Gauss_quadrature. */

class WignerChebyshev2nd1DQuadrature : public Base1DQuadrature {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default constructor.
  /*! Constructs a Gaussian-Chebyshev quadrature (of second type) of order \c order,
   * in the interval <c>[minDomainValue,maxDomainValue]</c>. This method scales the
   * abscissas (positions) of the quadrature from the interval [-1,1] to
   * <c>[minDomainValue,maxDomainValue]</c>.*/

  WignerChebyshev2nd1DQuadrature(double       minDomainValue,
                                        double       maxDomainValue,
                                        unsigned int order);
  //! Destructor.
  ~WignerChebyshev2nd1DQuadrature();
  //@}

};

}  // End namespace QUESO

#endif // UQ_1D_1D_QUADRATURE_H
