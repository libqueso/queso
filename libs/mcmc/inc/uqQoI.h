/* uq/libs/mcmc/inc/uqQoI.h
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

#ifndef __UQ_QOI_H__
#define __UQ_QOI_H__

#include <string>
#include <iostream>

class uqQoIClass
{
public:
  uqQoIClass(const std::string& name);
 ~uqQoIClass();

  std::string  name   () const;
  void         setName(const std::string& name);
  void         print  (std::ostream& os) const;

private:
  std::string  m_name;
};

std::ostream& operator<<(std::ostream& os, const uqQoIClass& qoi);

#endif // __UQ_QOI_H__
