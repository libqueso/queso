/* uq/libs/mcmc/inc/uqStateSpace.h
 *
 * Copyright (C) 2008 The PECOS Team, http://www.ices.utexas.edu/centers/pecos
 *
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

#ifndef __UQ_STATE_SPACE_H__
#define __UQ_STATE_SPACE_H__

#include <uqFinDimLinearSpace.h>

template <class V, class M>
class uqStateSpaceClass : public uqFinDimLinearSpaceClass<V,M>
{
public:
           uqStateSpaceClass();
           uqStateSpaceClass(const uqEnvironmentClass& env,
                             const char*               prefix);
  virtual ~uqStateSpaceClass();

          unsigned int            dim                       () const;
  virtual void                    print                     (std::ostream& os) const;

protected:
          void                    defineMyOptions           (po::options_description& optionsDesc) const;
          void                    getMyOptionValues         (po::options_description& optionsDesc);

  po::options_description* m_optionsDesc;
  std::string m_option_help;
  std::string m_option_dim;

  using uqFinDimLinearSpaceClass<V,M>::m_env;
  using uqFinDimLinearSpaceClass<V,M>::m_dim;
};

template <class V, class M>
uqStateSpaceClass<V,M>::uqStateSpaceClass()
  :
  uqFinDimLinearSpaceClass<V,M>()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqStateSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqStateSpaceClass<V,M>::uqStateSpaceClass(
  const uqEnvironmentClass& env,
  const char*               prefix)
  :
  uqFinDimLinearSpaceClass<V,M>(env,prefix),
  m_optionsDesc                (NULL)
{
  //std::cout << "Entering uqStateSpaceClass<V,M>::constructor()"
  //          << std::endl;

  m_option_help = uqFinDimLinearSpaceClass<V,M>::m_prefix + "stateSpace_help";
  m_option_dim  = uqFinDimLinearSpaceClass<V,M>::m_prefix + "stateSpace_dim";

  m_optionsDesc = new po::options_description("State space options");
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  // Now that 'm_dim' has been set, construct Trilinos map
  this->constructMap();

  if (m_env.rank() == 0) std::cout << "After getting values of options with prefix '" <<  uqFinDimLinearSpaceClass<V,M>::m_prefix
                                   << "', state of uqStateSpaceClass object is:"
                                   << "\n" << *this
                                   << std::endl;

  //std::cout << "Leaving uqStateSpaceClass<V,M>::constructor()"
  //          << std::endl;
}

template <class V, class M>
uqStateSpaceClass<V,M>::~uqStateSpaceClass()
{
  //std::cout << "Entering uqStateSpaceClass<V,M>::destructor()"
  //          << std::endl;

  if (m_optionsDesc) delete m_optionsDesc;

  //std::cout << "Leaving uqStateSpaceClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
void
uqStateSpaceClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    (m_option_help.c_str(),                                              "produce help message for UQ state space")
    (m_option_dim.c_str(),  po::value<unsigned int>()->default_value(0), "Space dimension"                        )
  ;

  return;
}

template <class V, class M>
void
uqStateSpaceClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_dim.c_str())) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    m_dim = tmpMap[m_option_dim.c_str()].as<unsigned int>();
  }

  return;
}

template <class V, class M>
unsigned int
uqStateSpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
void
uqStateSpaceClass<V,M>::print(std::ostream& os) const
{
  os <<  uqFinDimLinearSpaceClass<V,M>::m_prefix << "dim = " << m_dim
     << std::endl;

  return;
}

template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqStateSpaceClass<V,M>& space)
{
  space.print(os);

  return os;
}
#endif // __UQ_STATE_SPACE_H__

