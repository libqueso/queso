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

#ifndef QUESO_ALGORITHM_FACTORY_H
#define QUESO_ALGORITHM_FACTORY_H

#include <queso/Factory.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/VectorSpace.h>
#include <queso/Algorithm.h>

namespace QUESO
{

/**
 * AlgorithmFactory class defintion.  Clients subclassing this for their own
 * algorithm (aka Metropolis-Hastings acceptance ratio) should implement
 * build_algorithm.
 */
class AlgorithmFactory : public Factory<Algorithm<GslVector, GslMatrix> >
{
public:
  /**
   * Constructor. Takes the name to be mapped.
   */
  AlgorithmFactory(const std::string & name)
    : Factory<Algorithm<GslVector, GslMatrix> >(name)
  {}

  /**
   * Destructor. (Empty.)
   */
  virtual ~AlgorithmFactory() {}

  static void set_environment(const BaseEnvironment & env)
  {
    m_env = &env;
  }

  static void set_tk(const BaseTKGroup<GslVector, GslMatrix> & tk)
  {
    m_tk = &tk;
  }

protected:
  virtual SharedPtr<Algorithm<GslVector, GslMatrix> >::Type build_algorithm() = 0;

  static const BaseEnvironment * m_env;
  static const BaseTKGroup<GslVector, GslMatrix> * m_tk;

private:
  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual SharedPtr<Algorithm<GslVector, GslMatrix> >::Type create();
};

inline
SharedPtr<Algorithm<GslVector, GslMatrix> >::Type
AlgorithmFactory::create()
{
  queso_require_msg(m_env, "ERROR: must call set_environment() before building alg!");
  queso_require_msg(m_tk, "ERROR: must call set_tk() before building alg!");

  SharedPtr<Algorithm<GslVector, GslMatrix> >::Type new_alg = this->build_algorithm();

  queso_assert(new_alg);

  return new_alg;
}

/**
 * AlgorithmFactoryImp implementation of AlgorithmFactory.  Implements an
 * algorithm factory for the standard Metropolis-Hastings algorithm (aka
 * acceptance ratio).
 */
template <class DerivedAlgorithm>
class AlgorithmFactoryImp : public AlgorithmFactory
{
public:
  AlgorithmFactoryImp(const std::string & name)
    : AlgorithmFactory(name)
  {}

  virtual ~AlgorithmFactoryImp() {}

private:
  virtual SharedPtr<Algorithm<GslVector, GslMatrix> >::Type build_algorithm()
  {
    SharedPtr<Algorithm<GslVector, GslMatrix> >::Type new_alg;
    new_alg.reset(new DerivedAlgorithm(*(this->m_env), *(this->m_tk)));
    return new_alg;
  }

};

} // namespace QUESO

#endif // QUESO_ALGORITHM_FACTORY_H
