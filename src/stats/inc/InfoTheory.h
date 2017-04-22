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

#ifndef UQ_INFO_THEORY_H
#define UQ_INFO_THEORY_H

#include <queso/Defines.h>
#ifdef QUESO_HAS_ANN

#include <ANN/ANN.h>
#include <ANN/ANNx.h>

// TODO: create InfoTheoryOptions
#define UQ_INFTH_ANN_NO_SMP        10000
#define UQ_INFTH_ANN_EPS           0.0
#define UQ_INFTH_ANN_KNN           6

namespace QUESO {

/*!
 * For each point q in dataX, this performs a k-Nearest Neighbours search
 * of the points in dataY closest to q.  The array distsXY is filled with the
 * distance from q to the k-th nearest neighbour.
 */
void distANN_XY( const ANNpointArray dataX, const ANNpointArray dataY,
                 double* distsXY,
                 unsigned int dimX, unsigned int dimY,
                 unsigned int xN, unsigned int yN,
                 unsigned int k, double eps );

/*!
 * Normalises dataXY by its sample mean and covariance.
 * dataX and dataY are populated with the marginals of dataXY.
 */
void normalizeANN_XY( ANNpointArray dataXY, unsigned int dimXY,
                      ANNpointArray dataX, unsigned int dimX,
                      ANNpointArray dataY, unsigned int dimY,
                      unsigned int N );

/*!
 * Not defined anywhere!
 */
void whiteningANN_X_Y( ANNpointArray dataX1, ANNpointArray dataX2,
                       unsigned int dimX, unsigned int N1, unsigned int N2 );

/*!
 * Computes the mutual information
 */
double computeMI_ANN( ANNpointArray dataXY,
                      unsigned int dimX, unsigned int dimY,
                      unsigned int k, unsigned int N, double eps );

/*!
 * Function: estimateMI_ANN (using a joint)
 * (Mutual Information)
 *
 * Computes the mutual information using a queso joint RV
 */
template <template <class P_V, class P_M> class RV, class P_V, class P_M>
double estimateMI_ANN( const RV<P_V,P_M>& jointRV,
                       const unsigned int xDimSel[], unsigned int dimX,
                       const unsigned int yDimSel[], unsigned int dimY,
                       unsigned int k, unsigned int N, double eps );

/*!
 * Function: estimateMI_ANN (using two seperate RVs)
 * (Mutual Information)
 */
template <class P_V, class P_M,
  template <class P_V, class P_M> class RV_1,
  template <class P_V, class P_M> class RV_2>
double estimateMI_ANN( const RV_1<P_V,P_M>& xRV,
                       const RV_2<P_V,P_M>& yRV,
                       const unsigned int xDimSel[], unsigned int dimX,
                       const unsigned int yDimSel[], unsigned int dimY,
                       unsigned int k, unsigned int N, double eps );

/*!
 * Function: estimateKL_ANN
 * (Kullback-Leibler divergence)
 *
 * Computes the Kullback-Leibler divergence between two queso random variables
 * xRV and yRV.
 */
template <class P_V, class P_M,
  template <class P_V, class P_M> class RV_1,
  template <class P_V, class P_M> class RV_2>
double estimateKL_ANN( RV_1<P_V,P_M>& xRV,
                       RV_2<P_V,P_M>& yRV,
                       unsigned int xDimSel[], unsigned int dimX,
                       unsigned int yDimSel[], unsigned int dimY,
                       unsigned int xN, unsigned int yN,
                       unsigned int k, double eps );

/*!
 * Function: estimateCE_ANN
 * (Cross Entropy)
 *
 * Estimates the cross-entropy of two queso random variables xRV and yRV
 */
template <class P_V, class P_M,
  template <class P_V, class P_M> class RV_1,
  template <class P_V, class P_M> class RV_2>
double estimateCE_ANN( RV_1<P_V,P_M>& xRV,
                       RV_2<P_V,P_M>& yRV,
                       unsigned int xDimSel[], unsigned int dimX,
                       unsigned int yDimSel[], unsigned int dimY,
                       unsigned int xN, unsigned int yN,
                       unsigned int k, double eps );

}  // End namespace QUESO

#endif // QUESO_HAS_ANN

#endif // UQ_INFO_THEORY_H
