//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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
// 
// $Id: uqInfoTheory.h 14612 2011-01-02 12:39:58Z gabriel $
//
//--------------------------------------------------------------------------

#ifndef UQ_INFO_THEORY_H
#define UQ_INFO_THEORY_H

#include <queso/Defines.h>
#ifdef QUESO_HAS_ANN

#include <ANN/ANN.h>
#include <ANN/ANNx.h>
#include <gsl/gsl_sf_psi.h>

// TODO: create InfoTheoryOptions
#define UQ_INFTH_ANN_NO_SMP        10000
#define UQ_INFTH_ANN_EPS           0.0
#define UQ_INFTH_ANN_KNN           6

namespace QUESO {

void distANN_XY( const ANNpointArray dataX, const ANNpointArray dataY, 
		 double* distsXY, 
		 unsigned int dimX, unsigned int dimY, 
		 unsigned int xN, unsigned int yN, 
		 unsigned int k, double eps );

void normalizeANN_XY( ANNpointArray dataXY, unsigned int dimXY,
		      ANNpointArray dataX, unsigned int dimX,
		      ANNpointArray dataY, unsigned int dimY,
		      unsigned int N );

void whiteningANN_X_Y( ANNpointArray dataX1, ANNpointArray dataX2, 
		       unsigned int dimX, unsigned int N1, unsigned int N2 );

double computeMI_ANN( ANNpointArray dataXY,
		      unsigned int dimX, unsigned int dimY,
		      unsigned int k, unsigned int N, double eps );

//*****************************************************
// Function: estimateMI_ANN (using a joint)
// (Mutual Information)
//*****************************************************
template<template <class P_V, class P_M> class RV, class P_V, class P_M>
double estimateMI_ANN( const RV<P_V,P_M>& jointRV, 
		       const unsigned int xDimSel[], unsigned int dimX,
		       const unsigned int yDimSel[], unsigned int dimY,
		       unsigned int k, unsigned int N, double eps )
{
  ANNpointArray dataXY;
  double MI_est;

  unsigned int dimXY = dimX + dimY;

  // Allocate memory
  dataXY = annAllocPts(N,dimXY);

  // Copy samples in ANN data structure
  P_V smpRV( jointRV.imageSet().vectorSpace().zeroVector() );
  for( unsigned int i = 0; i < N; i++ ) {
    // get a sample from the distribution
    jointRV.realizer().realization( smpRV );

    // copy the vector values in the ANN data structure
    for( unsigned int j = 0; j < dimX; j++ ) {
      dataXY[ i ][ j ] = smpRV[ xDimSel[j] ];
    }
    for( unsigned int j = 0; j < dimY; j++ ) {
      dataXY[ i ][ dimX + j ] = smpRV[ yDimSel[j] ];
    }
    // annPrintPt( dataXY[i], dimXY, std::cout ); std::cout << std::endl;
  }

  MI_est = computeMI_ANN( dataXY,
			  dimX, dimY,
			  k, N, eps );

  // Deallocate memory
  annDeallocPts( dataXY );

  return MI_est;
}

//*****************************************************
// Function: estimateMI_ANN (using two seperate RVs)
// (Mutual Information)
//*****************************************************
template<class P_V, class P_M,
  template <class P_V, class P_M> class RV_1,
  template <class P_V, class P_M> class RV_2>
double estimateMI_ANN( const RV_1<P_V,P_M>& xRV, 
		       const RV_2<P_V,P_M>& yRV, 
		       const unsigned int xDimSel[], unsigned int dimX,
		       const unsigned int yDimSel[], unsigned int dimY,
		       unsigned int k, unsigned int N, double eps )
{
  ANNpointArray dataXY;
  double MI_est;

  unsigned int dimXY = dimX + dimY;

  // Allocate memory
  dataXY = annAllocPts(N,dimXY);

  // Copy samples in ANN data structure
  P_V smpRV_x( xRV.imageSet().vectorSpace().zeroVector() );
  P_V smpRV_y( yRV.imageSet().vectorSpace().zeroVector() );

  for( unsigned int i = 0; i < N; i++ ) {
    // get a sample from the distribution
    xRV.realizer().realization( smpRV_x );
    yRV.realizer().realization( smpRV_y );

    // copy the vector values in the ANN data structure
    for( unsigned int j = 0; j < dimX; j++ ) {
      dataXY[ i ][ j ] = smpRV_x[ xDimSel[j] ];
    }
    for( unsigned int j = 0; j < dimY; j++ ) {
      dataXY[ i ][ dimX + j ] = smpRV_y[ yDimSel[j] ];
    }
    // annPrintPt( dataXY[i], dimXY, std::cout ); std::cout << std::endl;
  }

  MI_est = computeMI_ANN( dataXY,
			  dimX, dimY,
			  k, N, eps );

  // Deallocate memory
  annDeallocPts( dataXY );

  return MI_est;
}

//*****************************************************
// Function: estimateKL_ANN
// (Kullback-Leibler divergence)
//*****************************************************
template <class P_V, class P_M, 
  template <class P_V, class P_M> class RV_1,
  template <class P_V, class P_M> class RV_2>
double estimateKL_ANN( RV_1<P_V,P_M>& xRV, 
		       RV_2<P_V,P_M>& yRV, 
		       unsigned int xDimSel[], unsigned int dimX,
		       unsigned int yDimSel[], unsigned int dimY,
		       unsigned int xN, unsigned int yN,
		       unsigned int k, double eps )
{
  ANNpointArray dataX;
  ANNpointArray dataY;
  double* distsX;
  double* distsXY;
  double KL_est;

  // sanity check
  if( dimX != dimY ) {
    std::cout << "Error-KL: the dimensions should agree" << std::endl;
    std::exit(1);
  }

  // Allocate memory
  dataX = annAllocPts( xN, dimX );
  dataY = annAllocPts( yN, dimY );
  distsX = new double[xN];
  distsXY = new double[xN];
  
  // Copy X samples in ANN data structure
  P_V xSmpRV( xRV.imageSet().vectorSpace().zeroVector() );
  for( unsigned int i = 0; i < xN; i++ ) {
    // get a sample from the distribution
    xRV.realizer().realization( xSmpRV );
    // copy the vector values in the ANN data structure
    for( unsigned int j = 0; j < dimX; j++ ) {
      dataX[ i ][ j ] = xSmpRV[ xDimSel[j] ];
    }
  }

  // Copy Y samples in ANN data structure
  P_V ySmpRV( yRV.imageSet().vectorSpace().zeroVector() );
  for( unsigned int i = 0; i < yN; i++ ) {
    // get a sample from the distribution
    yRV.realizer().realization( ySmpRV );
    // copy the vector values in the ANN data structure
    for( unsigned int j = 0; j < dimY; j++ ) {
      dataY[ i ][ j ] = ySmpRV[ yDimSel[j] ];
    }
  }

  // Get distance to knn for each point
  distANN_XY( dataX, dataX, distsX, dimX, dimX, xN, xN, k+1, eps ); // k+1 because the 1st nn is itself
  distANN_XY( dataX, dataY, distsXY, dimX, dimY, xN, yN, k, eps );
  
  // Compute KL-divergence estimate
  double sum_log_ratio = 0.0;
  for( unsigned int i = 0; i < xN; i++ ) 
    {
      sum_log_ratio += log( distsXY[i] / distsX[i] );
    }
  KL_est = (double)dimX/(double)xN * sum_log_ratio + log( (double)yN / ((double)xN-1.0 ) );

  // Deallocate memory
  annDeallocPts( dataX );
  annDeallocPts( dataY );
  delete [] distsX;
  delete [] distsXY;

  return KL_est;
}


//*****************************************************
// Function: estimateCE_ANN
// (Cross Entropy)
//*****************************************************
template <class P_V, class P_M, 
  template <class P_V, class P_M> class RV_1,
  template <class P_V, class P_M> class RV_2>
double estimateCE_ANN( RV_1<P_V,P_M>& xRV, 
		       RV_2<P_V,P_M>& yRV, 
		       unsigned int xDimSel[], unsigned int dimX,
		       unsigned int yDimSel[], unsigned int dimY,
		       unsigned int xN, unsigned int yN,
		       unsigned int k, double eps )
{
  ANNpointArray dataX;
  ANNpointArray dataY;
  double* distsXY;
  double CE_est;
  ANNkd_tree* kdTree;

  // sanity check
  if( dimX != dimY ) {
    std::cout << "Error-CE: the dimensions should agree" << std::endl;
    std::exit(1);
  }

  // Allocate memory
  dataX = annAllocPts( xN, dimX );
  dataY = annAllocPts( yN, dimY );
  distsXY = new double[xN];
  kdTree = new ANNkd_tree( dataY, yN, dimY );
  
  // Copy X samples in ANN data structure
  P_V xSmpRV( xRV.imageSet().vectorSpace().zeroVector() );
  for( unsigned int i = 0; i < xN; i++ ) {
    // get a sample from the distribution
    xRV.realizer().realization( xSmpRV );
    // copy the vector values in the ANN data structure
    for( unsigned int j = 0; j < dimX; j++ ) {
      dataX[ i ][ j ] = xSmpRV[ xDimSel[j] ];
    }
  }

  // Copy Y samples in ANN data structure
  P_V ySmpRV( yRV.imageSet().vectorSpace().zeroVector() );
  for( unsigned int i = 0; i < yN; i++ ) {
    // get a sample from the distribution
    yRV.realizer().realization( ySmpRV );
    // copy the vector values in the ANN data structure
    for( unsigned int j = 0; j < dimY; j++ ) {
      dataY[ i ][ j ] = ySmpRV[ yDimSel[j] ];
    }
  }

  // Get distance to knn for each point
  distANN_XY( dataX, dataY, distsXY, dimX, dimY, xN, yN, k, eps );
  kdTree = new ANNkd_tree( dataY, yN, dimY );

  // Compute cross entropy estimate
  double sum_log = 0.0;
  for( unsigned int i = 0; i < xN; i++ ) 
    {
      sum_log += log( distsXY[i] );
    }
  CE_est = (double)dimX/(double)xN * sum_log + log( (double)yN ) - gsl_sf_psi_int( k );

  // Deallocate memory
  annDeallocPts( dataX );
  annDeallocPts( dataY );
  delete [] distsXY;

  return CE_est;
}

}  // End namespace QUESO

#endif // QUESO_HAS_ANN

#endif // UQ_INFO_THEORY_H
