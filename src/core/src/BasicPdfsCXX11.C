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

#include <queso/BasicPdfsCXX11.h>

#ifdef QUESO_HAVE_CXX11

#include <cmath>

namespace QUESO {

BasicPdfsCXX11::BasicPdfsCXX11(int worldRank)
  :
  BasicPdfsBase(worldRank)
{
}

BasicPdfsCXX11::~BasicPdfsCXX11()
{
}

double
BasicPdfsCXX11::betaPdfActualValue(double x, double alpha, double beta) const
{
  // C++11 doesn't have a beta distribution, so we just implement one.

  double ln_gamma_apb;
  double ln_gamma_a;
  double ln_gamma_b;

  if (x < 0 || x > 1) {  // Outside the support; return zero.
    return 0.0;
  }
  else {
    if (x == 0.0 || x == 1.0) {  // On the boundary
      if (alpha > 1.0 && beta > 1.0) {  // Not blowing up
        return 0.0;
      }
      else {  // Might be blowing up, so let it.
        ln_gamma_apb = std::lgamma(alpha + beta);
        ln_gamma_a = std::lgamma(alpha);
        ln_gamma_b = std::lgamma(beta);

        return std::exp(ln_gamma_apb - ln_gamma_a - ln_gamma_b) *
          std::pow(x, alpha - 1) * std::pow(1 - x, beta - 1);
      }
    }
    else {  // Compute beta pdf in log space.
      ln_gamma_apb = std::lgamma(alpha + beta);
      ln_gamma_a = std::lgamma(alpha);
      ln_gamma_b = std::lgamma(beta);

      return std::exp(ln_gamma_apb - ln_gamma_a - ln_gamma_b +
          (alpha - 1) * std::log(x) + (beta - 1) * std::log1p(-x));
    }
  }
}

double
BasicPdfsCXX11::gammaPdfActualValue(double x, double a, double b) const
{
  if (x < 0) {  // Out of bounds.
    return 0;
  }
  else {
    if (x == 0.0) {  // On the boundary.
      if (a > 1.0) {  // Not blowing up.
        return 0.0;
      }
      else {
        double ln_gamma_a = std::lgamma(a);
        return std::exp(-ln_gamma_a - a * std::log(b)) * std::pow(x, a - 1) *
          std::exp(- x / b);
      }
    }
    else {
      double ln_gamma_a = std::lgamma(a);
      return std::exp(-ln_gamma_a - a * std::log(b) + (a - 1) * std::log(x) -
          (x / b));
    }
  }
}

}  // End namespace QUESO

#endif  // QUESO_HAVE_CXX11
