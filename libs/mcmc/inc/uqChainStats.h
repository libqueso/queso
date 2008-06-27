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

