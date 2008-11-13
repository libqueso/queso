/* uq/libs/queso/inc/uqTK.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_TRANSITION_KERNEL_H__
#define __UQ_TRANSITION_KERNEL_H__

#include <uqVectorSpace.h>
#include <uqTKPdf.h>
#include <uqTKRealizer.h>

//*****************************************************
// Base class
//*****************************************************
template<class V, class M>
class uqBaseTKClass {
public:
  uqBaseTKClass();
  uqBaseTKClass(const char*                    prefix,
                const uqVectorSpaceClass<V,M>& vectorSpace);
  virtual ~uqBaseTKClass();

  const   uqBaseEnvironmentClass&     env        () const;
  const   uqBaseTKPdfClass     <V,M>& pdf        () const;
  const   uqBaseTKRealizerClass<V,M>& realizer   () const;

  virtual void setPdf     (const uqBaseTKPdfClass     <V,M>& tKPdf     ) = 0;
  virtual void setRealizer(const uqBaseTKRealizerClass<V,M>& tKRealizer) = 0;

protected:
  const   uqEmptyEnvironmentClass*    m_emptyEnv;
  const   uqBaseEnvironmentClass&     m_env;
          std::string                 m_prefix;
  const   uqBaseTKPdfClass     <V,M>* m_tKPdf;
  const   uqBaseTKRealizerClass<V,M>* m_tKRealizer;
};

template<class V, class M>
uqBaseTKClass<V,M>::uqBaseTKClass()
  :
  m_emptyEnv  (new uqEmptyEnvironmentClass()),
  m_env       (*m_emptyEnv),
  m_prefix    (""),
  m_tKPdf     (NULL),
  m_tKRealizer(NULL)
{
}

template<class V, class M>
uqBaseTKClass<V,M>::uqBaseTKClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace)
  :
  m_emptyEnv  (NULL),
  m_env       (vectorSpace.env()),
  m_prefix    (prefix),
  m_tKPdf     (NULL),
  m_tKRealizer(NULL)
{
}

template<class V, class M>
uqBaseTKClass<V,M>::~uqBaseTKClass()
{
  if (m_emptyEnv) delete m_emptyEnv;
}

template<class V, class M>
const uqBaseEnvironmentClass&
uqBaseTKClass<V,M>::env() const
{
  return m_env;
}

template<class V, class M>
const uqBaseTKPdfClass<V,M>&
uqBaseTKClass<V,M>::pdf() const
{
  UQ_FATAL_TEST_MACRO(m_tKPdf == NULL,
                      m_env.rank(),
                      "uqBaseTKClass<V,M>::pdf()",
                      "m_tKPdf is NULL");

  return *m_tKPdf;
}

template<class V, class M>
const uqBaseTKRealizerClass<V,M>&
uqBaseTKClass<V,M>::realizer() const
{
  UQ_FATAL_TEST_MACRO(m_tKRealizer == NULL,
                      m_env.rank(),
                      (std::string)("uqBaseTKClass<V,M>::realizer(), prefix=")+m_prefix,
                      "m_tKRealizer is NULL");

  return *m_tKRealizer;
}

//*****************************************************
// Langevin class
//*****************************************************
template<class V, class M>
class uqLangevinTKClass : public uqBaseTKClass<V,M> {
public:
  uqLangevinTKClass();
  uqLangevinTKClass(const char*                    prefix,
                    const uqVectorSpaceClass<V,M>& vectorSpace);
 ~uqLangevinTKClass();

  void setPdf     (const uqBaseTKPdfClass     <V,M>& tKPdf     );
  void setRealizer(const uqBaseTKRealizerClass<V,M>& tKRealizer);

private:
  using uqBaseTKClass<V,M>::m_env;
  using uqBaseTKClass<V,M>::m_prefix;
  using uqBaseTKClass<V,M>::m_tKPdf;
  using uqBaseTKClass<V,M>::m_tKRealizer;
};

template<class V, class M>
uqLangevinTKClass<V,M>::uqLangevinTKClass()
  :
  uqBaseTKClass<V,M>()
{
}

template<class V, class M>
uqLangevinTKClass<V,M>::uqLangevinTKClass(
  const char*                    prefix,
  const uqVectorSpaceClass<V,M>& vectorSpace)
  :
  uqBaseTKClass<V,M>(prefix,vectorSpace)
{
  m_tKPdf = NULL; //new uqLangevinTKPdfClass<V,M>(m_prefix.c_str(),
                  //                              vectorSpace,
                  //                              covMatrix);
  m_tKRealizer = NULL; // FIX ME: complete code
}

template<class V, class M>
uqLangevinTKClass<V,M>::~uqLangevinTKClass()
{
}

template<class V, class M>
void
uqLangevinTKClass<V,M>::setPdf(const uqBaseTKPdfClass<V,M>& tKPdf)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqLangevinTKClass<V,M>::setPdf()",
                      "it does not make sense to call such routine for this class");
  return;
}

template<class V, class M>
void
uqLangevinTKClass<V,M>::setRealizer(const uqBaseTKRealizerClass<V,M>& tKRealizer)
{
  UQ_FATAL_TEST_MACRO(true,
                      m_env.rank(),
                      "uqLangevinTKClass<V,M>::setRealizer()",
                      "it does not make sense to call such routine for this class");
  return;
}
#endif // __UQ_TRANSITION_KERNEL_H__
