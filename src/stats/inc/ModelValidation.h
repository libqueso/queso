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

#ifndef UQ_MODEL_VALIDATION_H
#define UQ_MODEL_VALIDATION_H

#include <queso/ValidationCycle.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*! \file ModelValidation.h
 * \brief A templated class for model validation of the example validationPyramid.
 *
 * \class ModelValidation
 * \brief A templated class for model validation of the example validationPyramid.
 *
 * Its derived class exPhysics1Validation enables comparison between the calibration
 * and validate stages. */

template <class P_V = GslVector, class P_M = GslMatrix, class Q_V = GslVector, class Q_M = GslMatrix>
class ModelValidation
{
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor.
  ModelValidation(const BaseEnvironment& env,
                         const char*                   prefix);

  //! Virtual destructor.
  virtual ~ModelValidation();
  //@}

  //! @name Misc methods
  //@{
  //! Runs calibration, validation and comparison stages. See template specialization.
  virtual void run() = 0;

  //! Access to the environment variable (m_env).
  const BaseEnvironment&                  env  () const;

  //! Access to the cycle (m_cycle).
  const ValidationCycle<P_V,P_M,Q_V,Q_M>& cycle() const;
  //@}

protected:
  const BaseEnvironment& m_env;
        std::string             m_prefix;

  ValidationCycle<P_V,P_M,Q_V,Q_M>* m_cycle;
};

}  // End namespace QUESO

#endif // UQ_MODEL_VALIDATION_H
