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

#include <queso/MarkovChainPositionData.h>

namespace QUESO {

class GslVector;
class GslMatrix;
class BaseEnvironment;
template <class V, class M> class BaseTKGroup;

template <class V = GslVector, class M = GslMatrix>
class Algorithm
{
public:
  Algorithm(const BaseEnvironment & env, const BaseTKGroup<V, M> & tk);
  ~Algorithm();

  //! Calculates the finite dimensional Metropolis-Hastings acceptance ratio.
  /*!
   * tk_pos_x is the position of the tk when evaluating for x
   * tk_pos_y is the position of the tk when evaluating for y
   *
   * This method is called by the delayed rejection procedure.
   */
  double acceptance_ratio(MarkovChainPositionData<V> x,
                          MarkovChainPositionData<V> y,
                          const V & tk_pos_x,
                          const V & tk_pos_y);
private:
  const BaseEnvironment & m_env;
  const BaseTKGroup<V, M> & m_tk;
};

}  // End namespace QUESO
