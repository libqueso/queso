/* uq/libs/mcmc/src/uqQoI.C
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
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

#include <uqQoI.h>

uqQoIClass::uqQoIClass(
  const std::string& name)
  :
  m_name                (name)
{
}

uqQoIClass::~uqQoIClass()
{
}

std::string
uqQoIClass::name() const
{
  return m_name;
}

void
uqQoIClass::setName(const std::string& name)
{
  m_name = name;
  return;
}

void
uqQoIClass::print(std::ostream& os) const
{
  os << this->m_name;

  return;
}

std::ostream&
operator<<(std::ostream& os, const uqQoIClass& qoi)
{
  qoi.print(os);

  return os;
}
