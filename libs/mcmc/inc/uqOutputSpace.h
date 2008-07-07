/* uq/libs/mcmc/inc/uqOutputSpace.h
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

#ifndef __UQ_OUTPUT_SPACE_H__
#define __UQ_OUTPUT_SPACE_H__

#include <uqFinDimLinearSpace.h>

template <class V, class M>
class uqOutputSpaceClass : public uqFinDimLinearSpaceClass<V,M>
{
public:
           uqOutputSpaceClass();
           uqOutputSpaceClass(const uqEnvironmentClass& env);
  virtual ~uqOutputSpaceClass();

          unsigned int            dim                       () const;
  virtual void                    print                     (std::ostream& os) const;

protected:
          void                    defineMyOptions           (po::options_description& optionsDesc) const;
          void                    getMyOptionValues         (po::options_description& optionsDesc);

  po::options_description* m_optionsDesc;

  using uqFinDimLinearSpaceClass<V,M>::m_env;
  using uqFinDimLinearSpaceClass<V,M>::m_dim;
};

template <class V, class M>
uqOutputSpaceClass<V,M>::uqOutputSpaceClass()
  :
  uqFinDimLinearSpaceClass<V,M>()
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqOutputSpaceClass<V,M>::constructor(), default",
                      "should not be used by user");
}

template <class V, class M>
uqOutputSpaceClass<V,M>::uqOutputSpaceClass(
  const uqEnvironmentClass& env)
  :
  uqFinDimLinearSpaceClass<V,M>(env),
  m_optionsDesc                (NULL)
{
  //std::cout << "Entering uqOutputSpaceClass<V,M>::constructor()"
  //          << std::endl;

  m_optionsDesc = new po::options_description("Output space options");
  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "After getting option values, state of uqOutputSpaceClass object is:"
                                   << "\ndimension = " << m_dim
                                   << "\n"
                                   << std::endl;

  //std::cout << "Leaving uqOutputSpaceClass<V,M>::constructor()"
  //          << std::endl;
}

template <class V, class M>
uqOutputSpaceClass<V,M>::~uqOutputSpaceClass()
{
  //std::cout << "Entering uqOutputSpaceClass<V,M>::destructor()"
  //          << std::endl;

  if (m_optionsDesc) delete m_optionsDesc;

  //std::cout << "Leaving uqOutputSpaceClass<V,M>::destructor()"
  //          << std::endl;
}

template <class V, class M>
void
uqOutputSpaceClass<V,M>::defineMyOptions(
  po::options_description& optionsDesc) const
{
  m_optionsDesc->add_options()
    ("uqOutputSpace_help",                                              "produce help message for UQ PS")
    ("uqOutputSpace_dim",  po::value<unsigned int>()->default_value(0), "Space dimension"               )
  ;

  return;
}

template <class V, class M>
void
uqOutputSpaceClass<V,M>::getMyOptionValues(po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count("uqOutputSpace_help")) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count("uqOutputSpace_dim")) {
    const po::variables_map& tmpMap = m_env.allOptionsMap();
    m_dim = tmpMap["uqOutputSpace_dim"].as<unsigned int>();
  }

  return;
}

template <class V, class M>
unsigned int
uqOutputSpaceClass<V,M>::dim() const
{
  return m_dim;
}

template <class V, class M>
void
uqOutputSpaceClass<V,M>::print(std::ostream& os) const
{
  os << "m_dim = " << m_dim
     << std::endl;

  return;
}

template<class V, class M>
std::ostream&
operator<<(std::ostream& os, const uqOutputSpaceClass<V,M>& space)
{
  space.print(os);

  return os;
}
#endif // __UQ_OUTPUT_SPACE_H__

