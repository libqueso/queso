/* uq/libs/mcmc/inc/uqChainStats.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos

 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_CHAIN_STATS_H__
#define __UQ_CHAIN_STATS_H__

#include <iostream>

template <class V>
void
uqChainStats(
  const std::vector<const V*>& chain,
  V&                           mean,
  V&                           std)
{
  bool bRC = ((chain[0]->size() == mean.size()) &&
              (chain[0]->size() == std.size() ));
  UQ_FATAL_TEST_MACRO(bRC == false,
                      chain[0]->env().rank(),
                      "uqChainStats<V>()",
                      "incompatible sizes on input vectors");

  double doubleChainSize = (double) chain.size();

  mean.cwSet(0.);
  for (unsigned int i = 0; i < chain.size(); ++i) {
    for (unsigned int j = 0; j < mean.size(); ++j) {
      mean[j] += (*chain[i])[j]/doubleChainSize;
    }
  }

  std.cwSet(0.);
  for (unsigned int i = 0; i < chain.size(); ++i) {
    for (unsigned int j = 0; j < std.size(); ++j) {
      double diff = (*chain[i])[j] - mean[j];
      std[j] += diff*diff;
    }
  }
  for (unsigned int j = 0; j < std.size(); ++j) {
    std[j] = sqrt(std[j])/(doubleChainSize - 1.);
  }

  return;
}

#endif // __UQ_CHAIN_STATS_H__

