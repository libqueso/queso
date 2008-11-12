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

  virtual void setTKPdf     (const uqBaseTKPdfClass     <V,M>& tKPdf     ) = 0;
  virtual void setTKRealizer(const uqBaseTKRealizerClass<V,M>& tKRealizer) = 0;

protected:
  const   uqEmptyEnvironmentClass     m_emptyEnv;
  const   uqEnvironmentClass&         m_env;
          std::string                 m_prefix;
  const   uqBaseTKPdfClass     <V,M>* m_tKPdf;
  const   uqBaseTKRealizerClass<V,M>* m_tKRealizer;
};

template<class V, class M>
uqBaseTKClass<V,M>::uqBaseTKClass()
  :
  m_env       (m_emptyEnv),
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
  m_env       (vectorSpace.env()),
  m_prefix    (prefix),
  m_tKPdf     (NULL),
  m_tKRealizer(NULL)
{
}

template<class V, class M>
uqBaseTKClass<V,M>::~uqBaseTKClass()
{
}

//template<class V, class M>
//void
//uqBaseTKClass<V,M>::setTKPdf(uqBaseTKPdfClass<V,M>& tKPdf)
//{
//  return;
//}
#endif // __UQ_TRANSITION_KERNEL_H__
