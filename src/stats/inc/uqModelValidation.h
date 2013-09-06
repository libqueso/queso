//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_MODEL_VALIDATION_H__
#define __UQ_MODEL_VALIDATION_H__

#include <uqValidationCycle.h>

namespace QUESO {

/*! \file uqModelValidation.h
 * \brief A templated class for model validation of the example validationPyramid.
 *
 * \class uqModelValidationClass
 * \brief A templated class for model validation of the example validationPyramid. 
 * 
 * Its derived class exPhysics1ValidationClass enables comparison between the calibration
 * and validate stages. */

template <class P_V,class P_M,class Q_V,class Q_M>
class uqModelValidationClass
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  uqModelValidationClass(const uqBaseEnvironmentClass& env,
                         const char*                   prefix);
  
  //! Virtual destructor.
  virtual ~uqModelValidationClass();
  //@}
  
  //! @name Misc methods
  //@{
  //! Runs calibration, validation and comparison stages. See template specialization.
  virtual void run() = 0;

  //! Access to the environment variable (m_env).
  const uqBaseEnvironmentClass&                  env  () const;
  
  //! Access to the cycle (m_cycle).
  const uqValidationCycleClass<P_V,P_M,Q_V,Q_M>& cycle() const;
  //@}
  
protected:
  const uqBaseEnvironmentClass& m_env;
        std::string             m_prefix;

  uqValidationCycleClass<P_V,P_M,Q_V,Q_M>* m_cycle;
};
// Constructor -------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
uqModelValidationClass<P_V,P_M,Q_V,Q_M>::uqModelValidationClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  m_env   (env),
  m_prefix((std::string)(prefix) + ""),
  m_cycle (NULL)
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqModelValidationClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqModelValidationClass<P_V,P_M,Q_V,Q_M>::constructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  return;
}
// Destructor---------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
uqModelValidationClass<P_V,P_M,Q_V,Q_M>::~uqModelValidationClass()
{
  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Entering uqModeValidation::destructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }

  if (m_cycle) delete m_cycle;

  if (m_env.subDisplayFile()) {
    *m_env.subDisplayFile() << "Leaving uqModeValidation::destructor()"
                           << ": prefix = " << m_prefix
                           << std::endl;
  }
}
// Misc methods-------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const uqBaseEnvironmentClass&
uqModelValidationClass<P_V,P_M,Q_V,Q_M>::env() const
{
  return m_env;
}
//--------------------------------------------------
template <class P_V,class P_M,class Q_V,class Q_M>
const uqValidationCycleClass<P_V,P_M,Q_V,Q_M>&
uqModelValidationClass<P_V,P_M,Q_V,Q_M>::cycle() const
{
  return *m_cycle;
}

}  // End namespace QUESO

#endif // __UQ_MODEL_VALIDATION_H__
