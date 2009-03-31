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

#ifndef __UQ_MISCELLANEOUS_H__
#define __UQ_MISCELLANEOUS_H__

#include <gsl/gsl_rng.h>
#include <string>
#include <vector>

void   uqMiscReadDoublesFromString      (const std::string&   inputString,
                                         std::vector<double>& outputDoubles);
void   uqMiscReadWordsFromString        (const std::string&        inputString,
                                         std::vector<std::string>& outputWords);
//void   uqMiscExtractDoubleFromString    (std::string& inputString,
//                                         double&      outputDouble);
//void   uqMiscExtractWordFromString      (std::string& inputString,
//                                         std::string& outputWord);
int    uqMiscReadStringAndDoubleFromFile(std::ifstream& ifs,
                                         std::string&   termString,
                                         double*        termValue);
int    uqMiscReadCharsAndDoubleFromFile (std::ifstream& ifs,
                                         std::string&   termString,
                                         double*        termValue,
                                         bool&          endOfLineAchieved);
double uqMiscGammar                     (double         a,
                                         double         b,
                                         const gsl_rng* rng);
double uqMiscGetEllapsedSeconds         (struct timeval *timeval0);

double uqMiscHammingWindow              (unsigned int N, unsigned int j);

double uqMiscGaussianDensity            (double x, double mu, double sigma);

template <class V>
void
uqMiscComputePositionsBetweenMinMax(V                minValues,
                                    V                maxValues,
                                    std::vector<V*>& positions)
{
  double factor = 0.5;
  switch (positions.size()) {
    case 0:
      // Do nothing
    break;

    case 1:
      positions[0] = new V((1. - factor) * minValues + factor * maxValues);
    break;

    default:
      for (unsigned int i = 0; i < positions.size(); ++i) {
        factor = ((double) i)/(((double) positions.size()) - 1.);
        positions[i] = new V((1. - factor) * minValues + factor * maxValues);
      }
    break;
  }

  return;
}

#endif // __UQ_MISCELLANEOUS_H__
