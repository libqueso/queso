/* libs/basic/src/uqDefines.C
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

#include <uqDefines.h>

int uqMyRank() {
  int result = 0;
#ifdef __APPL_USES_TRILINOS__
  int iRC;
  iRC = MPI_Comm_rank(MPI_COMM_WORLD,&result);
#endif
  return result;
}
