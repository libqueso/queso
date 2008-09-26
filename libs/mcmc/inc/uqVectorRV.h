/* uq/libs/mcmc/inc/uqVectorRV.h
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

#ifndef __UQ_VECTOR_RV_H__
#define __UQ_VECTOR_RV_H__

#include <uqBasicScalarRV.h>
#include <uqProbDensity.h>
#include <uqRealizer.h>
#include <uqVectorSpace.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqVectorRVClass {
public:
  uqVectorRVClass(const uqEnvironmentClass&           env,
                  const char*                         prefix,
                  const uqVectorSpaceClass     <V,M>& imageSpace,
                  const uqProbDensity_BaseClass<V,M>* probDensity,
                  const uqRealizer_BaseClass   <V,M>* realizer);
  virtual ~uqVectorRVClass();

  const   uqProbDensity_BaseClass<V,M>& probDensity               ()       const;
  const   uqRealizer_BaseClass   <V,M>& realizer                  ()       const;
          void                          realization               (V& vec) const;

          int                           setParameter              (unsigned int       paramId,
                                                                   const std::string& name         = ".",
                                                                   double             initialValue = 0.,
                                                                   double             minValue     = -INFINITY,
                                                                   double             maxValue     = INFINITY,
                                                                   double             priorMu      = 0.,
                                                                   double             priorSigma   = INFINITY);
  const   uqBasicScalarRVClass&         parameter                 (unsigned int paramId) const;
  const   V&                            initialValues             () const;
  const   V&                            minValues                 () const;
  const   V&                            maxValues                 () const;
  const   V&                            priorMuValues             () const;
  const   V&                            priorSigmaValues          () const;
  const   std::vector<std::string>&     componentsNames           () const;

          bool                          outOfBounds               (const V& v) const;

  virtual void                          print                     (std::ostream& os) const;
          void                          printParameterNames       (std::ostream& os, bool printHorizontally) const; // See template specialization

protected:
          void                          defineMyOptions           (po::options_description& optionsDesc) const;
          void                          getMyOptionValues         (po::options_description& optionsDesc);
          void                          readParametersFromSpecFile(std::string& specFileName);
          void                          resetValues               ();
          void                          createInitialValues       () const; // See template specialization
          void                          createMinValues           () const; // See template specialization
          void                          createMaxValues           () const; // See template specialization
          void                          createPriorMuValues       () const; // See template specialization
          void                          createPriorSigmaValues    () const; // See template specialization
          void                          createComponentsNames     () const; // See template specialization

  const   uqEnvironmentClass&                m_env;
          std::string                        m_prefix;
  const   uqVectorSpaceClass     <V,M>*      m_imageSpace;
  const   uqProbDensity_BaseClass<V,M>*      m_probDensity;
  const   uqRealizer_BaseClass   <V,M>*      m_realizer;

          po::options_description*           m_optionsDesc;
          std::string                        m_option_help;
          std::string                        m_option_dim;
          std::string                        m_option_specificationFile;

          std::vector<uqBasicScalarRVClass*> m_parameters; // FIXME: will need to be a parallel vector in case of a very large number of parameters
          uqBasicScalarRVClass               m_dummyParameter;
          mutable V*                         m_initialValues;
          mutable V*                         m_minValues;
          mutable V*                         m_maxValues;
          mutable V*                         m_priorMuValues;
          mutable V*                         m_priorSigmaValues;
          mutable std::vector<std::string>   m_componentsNames;
};

template<class V, class M>
uqVectorRVClass<V,M>::uqVectorRVClass(
  const uqEnvironmentClass&           env,
  const char*                         prefix,
  const uqVectorSpaceClass     <V,M>& imageSpace,
  const uqProbDensity_BaseClass<V,M>* probDensity,
  const uqRealizer_BaseClass   <V,M>* realizer)
  :
  m_env                     (env),
  m_prefix                  (prefix),
  m_imageSpace              (&imageSpace),
  m_probDensity             (probDensity),
  m_realizer                (realizer),
  m_optionsDesc             (new po::options_description("Parameter space options")),
  m_option_help             (m_prefix + "help"),
  m_option_dim              (m_prefix + "dim"),
  m_option_specificationFile(m_prefix + "specificationFile"),
  m_parameters              (0),//,NULL),
  m_dummyParameter          (".",0.),
  m_initialValues           (NULL),
  m_minValues               (NULL),
  m_maxValues               (NULL),
  m_priorMuValues           (NULL),
  m_priorSigmaValues        (NULL),
  m_componentsNames         (0)
{
}

template<class V, class M>
uqVectorRVClass<V,M>::~uqVectorRVClass()
{
}

template<class V, class M>
const uqProbDensity_BaseClass<V,M>&
uqVectorRVClass<V,M>::probDensity() const
{
  UQ_FATAL_TEST_MACRO(m_probDensity == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqVectorRVClass<V,M>::probDensity()",
                      "m_probDensity is NULL");

  return *m_probDensity;
}

template<class V, class M>
const uqRealizer_BaseClass<V,M>&
uqVectorRVClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqVectorRVClass<V,M>::realizer()",
                      "m_realizer is NULL");

  return *m_realizer;
}

template<class V, class M>
void
uqVectorRVClass<V,M>::realization(V& vec) const
{
  UQ_FATAL_TEST_MACRO(m_realizer == NULL,
                      UQ_UNAVAILABLE_RANK,
                      "uqVectorRVClass<V,M>::realization()",
                      "m_realizer is NULL");

  m_realizer->nextSample(vec);

  return;
}

#endif // __UQ_VECTOR_RV_H__
