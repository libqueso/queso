/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008,2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Basic API: Internal header
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

using namespace std;

namespace QUESO_Basic_API 
{

#ifdef _QUESO_Basic_API_DEF
#define GLOBAL
#else
#define GLOBAL extern
#endif

#ifdef _QUESO_Basic_API_DEF
  QUESO_Basic_Class *_QUESO_Basic = NULL; 
#else
  extern QUESO_Basic_Class *_QUESO_Basic;
#endif



}
