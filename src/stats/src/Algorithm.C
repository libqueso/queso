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

#include <queso/Environment.h>
#include <queso/math_macros.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/Algorithm.h>
#include <queso/TKGroup.h>
#include <queso/InvLogitGaussianJointPdf.h>

namespace QUESO {

template <class V, class M>
Algorithm<V, M>::Algorithm(const BaseEnvironment & env,
    const BaseTKGroup<V, M> & tk)
  :
    m_env(env),
    m_tk(tk)
{
}

template <class V, class M>
Algorithm<V, M>::~Algorithm()
{
}

template <class V, class M>
double
Algorithm<V, M>::acceptance_ratio(
    MarkovChainPositionData<V> x,
    MarkovChainPositionData<V> y,
    const V & tk_pos_x,
    const V & tk_pos_y)
{
  double alphaQuotient = 0.;
  if ((x.outOfTargetSupport() == false) &&
      (y.outOfTargetSupport() == false)) {
    if ((x.logTarget() == -INFINITY) ||
        (x.logTarget() ==  INFINITY) ||
        ( queso_isnan(x.logTarget())      )) {
      std::cerr << "WARNING In Algorithm<V,M>::alpha(x,y)"
                << ", worldRank "       << m_env.worldRank()
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ": x.logTarget() = " << x.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else if ((y.logTarget() == -INFINITY           ) ||
             (y.logTarget() ==  INFINITY           ) ||
             ( queso_isnan(y.logTarget()) )) {
      std::cerr << "WARNING In Algorithm<V,M>::alpha(x,y)"
                << ", worldRank "       << m_env.worldRank()
                << ", fullRank "        << m_env.fullRank()
                << ", subEnvironment "  << m_env.subId()
                << ", subRank "         << m_env.subRank()
                << ", inter0Rank "      << m_env.inter0Rank()
                << ": y.logTarget() = " << y.logTarget()
                << ", x.values() = "    << x.vecValues()
                << ", y.values() = "    << y.vecValues()
                << std::endl;
    }
    else {
      double yLogTargetToUse = y.logTarget();

      if (m_tk.symmetric()) {
        alphaQuotient = std::exp(yLogTargetToUse - x.logTarget());

        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 3            )) {
          *m_env.subDisplayFile() << "In Algorithm<V,M>::alpha(x,y)"
                                 << ": symmetric proposal case"
                                 << ", x = "               << x.vecValues()
                                 << ", y = "               << y.vecValues()
                                 << ", yLogTargetToUse = " << yLogTargetToUse
                                 << ", x.logTarget() = "   << x.logTarget()
                                 << ", alpha = "           << alphaQuotient
                                 << std::endl;
        }
      }
      else {
        double qyx = m_tk.rv(tk_pos_x).pdf().lnValue(x.vecValues());
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 10           )) {
          *m_env.subDisplayFile() << m_tk.rv(tk_pos_x).pdf() << std::endl;
        }
        double qxy = m_tk.rv(tk_pos_y).pdf().lnValue(y.vecValues());
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 10           )) {
          *m_env.subDisplayFile() << m_tk.rv(tk_pos_y).pdf() << std::endl;
        }
        alphaQuotient = std::exp(yLogTargetToUse +
                                 qyx -
                                 x.logTarget() -
                                 qxy);
        if ((m_env.subDisplayFile()                   ) &&
            (m_env.displayVerbosity() >= 3            )) {
          *m_env.subDisplayFile() << "In Algorithm<V,M>::alpha(x,y)"
                                 << ": asymmetric proposal case"
                                 << ", x = "               << x.vecValues()
                                 << ", y = "               << y.vecValues()
                                 << ", yLogTargetToUse = " << yLogTargetToUse
                                 << ", q(y,x) = "          << qyx
                                 << ", x.logTarget() = "   << x.logTarget()
                                 << ", q(x,y) = "          << qxy
                                 << ", alpha = "           << alphaQuotient
                                 << std::endl;
        }
      }
    } // protection logic against logTarget values
  }
  else {
    if ((m_env.subDisplayFile()                   ) &&
        (m_env.displayVerbosity() >= 10           )) {
      *m_env.subDisplayFile() << "In Algorithm<V,M>::alpha(x,y)"
                             << ": x.outOfTargetSupport = " << x.outOfTargetSupport()
                             << ", y.outOfTargetSupport = " << y.outOfTargetSupport()
                             << std::endl;
    }
  }

  return std::min(1.,alphaQuotient);
}

template class Algorithm<GslVector, GslMatrix>;

}  // End namespace QUESO
