/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __EX_HYST_MODEL_H__
#define __EX_HYST_MODEL_H__

#include <queso/SequenceOfVectors.h>
#include <queso/GslMatrix.h>

void hystereticModel(
  const QUESO::BaseEnvironment& env,
  const QUESO::GslVector&       massInputVec,
  const QUESO::GslVector&       kInputVec,
  const QUESO::GslVector&       rInputVec,
  const QUESO::GslVector&       uInputVec,
  double                        rho,
  double                        gamma,
  const std::vector<double>&    a,
  std::vector<double>&          t,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& u,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& ud,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& udd,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& resfor,
  QUESO::SequenceOfVectors<QUESO::GslVector,QUESO::GslMatrix>& ru);

#endif
