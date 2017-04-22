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

#ifndef QUESO_FACTORY_H
#define QUESO_FACTORY_H

#include <queso/asserts.h>
#include <queso/SharedPtr.h>

#include <cstddef>
#include <map>
#include <string>
#include <iostream>

namespace QUESO
{

/**
 * Factory class defintion.
 */
template <class Base>
class Factory
{
protected:
  /**
   * Constructor. Takes the name to be mapped.
   */
  Factory(const std::string & name);

public:
  /**
   * Destructor. (Empty.)
   */
  virtual ~Factory() {}

  /**
   * Builds an object of type Base identified by name.
   */
  static typename SharedPtr<Base>::Type build(const std::string & name);

  /**
   * Create a Base class.  Force this to be implemented
   * later.
   */
  virtual typename SharedPtr<Base>::Type create() = 0;

protected:
  /**
   * Map from a name to a Factory<Base> * pointer.
   */
  static std::map<std::string, Factory<Base> *> & factory_map();
};

/**
 * Factory implementation class.
 */
template <class Derived, class Base>
class FactoryImp: public Factory<Base>
{
public:
  /**
   * Constructor.  Takes a name as input.
   */
  FactoryImp(const std::string & name) : Factory<Base>(name) { }

  /**
   * Destructor.  Empty.
   */
  ~FactoryImp() {}

private:
  /**
   * @returns a new object of type Derived.
   */
  virtual typename SharedPtr<Base>::Type create();
};

// -----------------------------------------------------
// Factory members
template <class Base>
inline
Factory<Base>::Factory(const std::string & name)
{
  // Make sure we haven't already added this name
  // to the map
  queso_assert(!factory_map().count(name));

  factory_map()[name] = this;
}

template <class Base>
inline
typename SharedPtr<Base>::Type Factory<Base>::build(const std::string & name)
{
  // name not found in the map
  if (!factory_map().count(name))
    {
      std::cerr << "Tried to build an unknown type: " << name << std::endl;

      std::cerr << "valid options are:" << std::endl;

      for (typename std::map<std::string,Factory<Base> *>::const_iterator
             it = factory_map().begin(); it != factory_map().end(); ++it)
        std::cerr << "  " << it->first << std::endl;

      queso_error_msg("Exiting...");

      // We'll never get here
      return typename SharedPtr<Base>::Type();
    }

  Factory<Base> * f = factory_map()[name];
  return typename SharedPtr<Base>::Type(f->create());
}

// Note - this cannot be inlined!
// template <class Base>
// std::map<std::string, Factory<Base> *> & Factory<Base>::factory_map()
// {
//   static std::map<std::string, Factory<Base> *> _factory_map;

//   return _factory_map;
// }

template <class Derived, class Base>
inline
typename SharedPtr<Base>::Type FactoryImp<Derived,Base>::create()
{
  return typename SharedPtr<Base>::Type(new Derived);
}

} // namespace QUESO

#endif // QUESO_FACTORY_H
