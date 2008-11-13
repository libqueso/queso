/* uq/examples/queso/inc/uqApplRoutines.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
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

#ifndef __UQ_APPL_ROUTINES_H__
#define __UQ_APPL_ROUTINES_H__

#include <uqEnvironment.h>

//*****************************************************
// User must provide an application routine, to be called by main().
//*****************************************************
template<class V, class M>
void
uqAppl(const uqBaseEnvironmentClass& env);

#endif // __UQ_APPL_ROUTINES_H__
