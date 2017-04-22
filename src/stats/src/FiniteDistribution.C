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

#include <queso/FiniteDistribution.h>
#include <queso/RngBase.h>

namespace QUESO {

// Default constructor -----------------------------
FiniteDistribution::FiniteDistribution(
  const BaseEnvironment& env,
  const char*                   prefix,
  const std::vector<double>&    inpWeights)
  :
  m_env    (env),
  m_prefix ((std::string)(prefix)+"fd_"),
  m_weights(inpWeights.size(),0.)
{
  queso_deprecated();

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Entering FiniteDistribution::constructor()"
                            << ": prefix = " << m_prefix
                            << ", inpWeights.size() = " << inpWeights.size()
                            << std::endl;
  }

  unsigned int numOfZeroWeights = 0;
  unsigned int numRareCases = 0;
  double sumCheck = 0.;
  unsigned int j = 0;
  m_map.empty(); // prudenci 2010-08-11
  for (unsigned int i = 0; i < inpWeights.size(); ++i) {
    double previousSum = sumCheck;
    sumCheck += inpWeights[i];
    if (sumCheck == previousSum) {
      numOfZeroWeights++;
    }
    else {
      if ((sumCheck - 1) > 1.e-8) {
        std::cerr << "In FiniteDistribution::constructor()"
                  << ": sumCheck - 1 = " << sumCheck - 1.
                  << std::endl;
      }
      queso_require_less_equal_msg((sumCheck - 1), 1.e-8, "weights sum is too bigger than 1.");

      if (sumCheck > 1.) sumCheck = 1.;
      m_weights[j] = inpWeights[i];
      std::pair<std::map<double,unsigned int>::iterator,bool> ret;
      ret = m_map.insert(std::map<double,unsigned int>::value_type(sumCheck,i));
      if (ret.second == true) {
        j++;
      }
      else {
        numRareCases++;
        if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
           *m_env.subDisplayFile() << "In FiniteDistribution::constructor()"
                                   << ": WARNING, map insertion failed"
                                   << std::endl;
        }
      }
    }
  }
  m_weights.resize(j,0.);

  if ((1 - sumCheck) > 1.e-8) {
    std::cerr << "In FiniteDistribution::constructor()"
              << ": 1 - sumCheck = " << 1. - sumCheck
              << std::endl;
  }
  queso_require_less_equal_msg((1 - sumCheck), 1.e-8, "weights sum is too smaller than 1.");


  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 3)) {
    *m_env.subDisplayFile() << "In FiniteDistribution::constructor()"
                            << ": inpWeights.size() = " << inpWeights.size()
                            << ", numOfZeroWeights = "  << numOfZeroWeights
                            << ", numRareCases = "      << numRareCases
                            << ", m_map.size() = "      << m_map.size()
                            << ", m_weights.size() = "  << m_weights.size()
                            << std::endl;
  }

  queso_require_equal_to_msg(inpWeights.size(), (m_weights.size()+numOfZeroWeights+numRareCases), "number of input weights was not conserved");

  queso_require_equal_to_msg(m_map.size(), m_weights.size(), "map and inpWeights have different sizes");

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 5)) {
    *m_env.subDisplayFile() << "Leaving FiniteDistribution::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
FiniteDistribution::~FiniteDistribution()
{
  queso_deprecated();

  m_map.empty();
  m_weights.clear();
}
// Misc methods--------------------------------------
const BaseEnvironment&
FiniteDistribution::env() const
{
  queso_deprecated();

  return m_env;
}
// Stats methods-------------------------------------
const std::vector<double>&
FiniteDistribution::weights() const
{
  queso_deprecated();

  return m_weights;
}
//---------------------------------------------------
unsigned int
FiniteDistribution::sample() const
{
  queso_deprecated();

  unsigned int result = 0;

  double aux = m_env.rngObject()->uniformSample();
  queso_require_msg(!((aux < 0) || (aux > 1.)), "invalid uniform");

  if (aux == 0.) {
    result = 0;
  }
  else if (aux == 1.) {
    result = m_map.find(aux)->second;
  }
  else {
    result = m_map.upper_bound(aux)->second;
    //if (m_map.upper_bound(aux)->second == 0) {
    //  result = 0;
    //}
    //else {
    //  result = m_map.upper_bound(aux)->second-1;
    //}
  }
#if 0 // WE insert 'i' in map, not 'j'. So, the tests below don't make sense
  if (result >= m_map.size()) {
    std::cerr << "In FiniteDistribution::sample()"
              << ": aux = "          << aux
              << ", m_map.size() = " << m_map.size()
              << ", result = "       << result
              << std::endl;
  }
  queso_require_less_msg(result, m_map.size(), "invalid result");
#endif

  return result;
}

}  // End namespace QUESO
