//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011 The PECOS Development Team
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

#ifndef __UQ_INFO_THEORY_H__
#define __UQ_INFO_THEORY_H__

#include <uqDefines.h>
#ifdef QUESO_HAS_ANN

#include <ANN/ANN.h>
#include <ANN/ANNx.h>
#include <gsl/gsl_sf_psi.h>

// TODO: create uqInfoTheoryOptions
#define UQ_INFTH_ANN_NO_SMP        10000
#define UQ_INFTH_ANN_EPS           0.0
#define UQ_INFTH_ANN_KNN           6

void distANN_XY( const ANNpointArray dataX, const ANNpointArray dataY, 
		 double* distsXY, 
		 unsigned int dimX, unsigned int dimY, 
		 unsigned int xN, unsigned int yN, 
		 unsigned int k, double eps );

void normalizeANN_XY( ANNpointArray dataXY, unsigned int dimXY,
		      ANNpointArray dataX, unsigned int dimX,
		      ANNpointArray dataY, unsigned int dimY,
		      unsigned int N );

//*****************************************************
// Function: estimateMI_ANN
//*****************************************************
template<template <class P_V, class P_M> class RV, class P_V, class P_M>
double estimateMI_ANN( const RV<P_V,P_M>& jointRV, 
		       const unsigned int xDimSel[], unsigned int dimX,
		       const unsigned int yDimSel[], unsigned int dimY,
		       unsigned int k, unsigned int N, double eps )
{
  ANNpointArray dataXY;
  ANNpointArray dataX, dataY;
  double* distsXY;
  double MI_est;
  ANNkd_tree* kdTreeX;
  ANNkd_tree* kdTreeY;

  unsigned int dimXY = dimX + dimY;

  // Allocate memory
  dataXY = annAllocPts(N,dimXY);
  dataX = annAllocPts(N,dimX);
  dataY = annAllocPts(N,dimY);
  distsXY = new double[N];

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

  // Normalize data and populate the marginals dataX, dataY
  normalizeANN_XY( dataXY, dimXY, dataX, dimX, dataY, dimY, N);

  // Get distance to knn for each point
  kdTreeX = new ANNkd_tree( dataX, N, dimX );
  kdTreeY = new ANNkd_tree( dataY, N, dimY );
  distANN_XY( dataXY, dataXY, distsXY, dimXY, dimXY, N, N, k, eps );

  // Compute mutual information
  double marginal_contrib = 0.0;
  for( unsigned int i = 0; i < N; i++ ) {
    // get the number of points within a specified radius
    int no_pts_X = kdTreeX->annkFRSearch( dataX[ i ], distsXY[ i ], 0, NULL, NULL, eps);
    int no_pts_Y = kdTreeY->annkFRSearch( dataY[ i ], distsXY[ i ], 0, NULL, NULL, eps);
    // digamma evaluations
    marginal_contrib += gsl_sf_psi_int( no_pts_X+1 ) + gsl_sf_psi_int( no_pts_Y+1 );
  }
  MI_est = gsl_sf_psi_int( k ) + gsl_sf_psi_int( N ) - marginal_contrib / (double)N;

  // Deallocate memory
  delete kdTreeX;
  delete kdTreeY;
  delete [] distsXY;
  annDeallocPts( dataX );  
  annDeallocPts( dataY );
  annDeallocPts( dataXY );

  return MI_est;
}


//*****************************************************
// Function: estimateKL_ANN
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

  // FIXME: here we have to make sure that  
  // the return value is a finite number
  // unsigned int xN = xRV->realizer().subPeriod();
  // unsigned int yN = yRV->realizer().subPeriod();

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
  for( unsigned int i = 0; i < xN; i++ ) {

    // FIXME: check if this can fail by given very large numbers
    double tmp_log = log( distsXY[i] / distsX[i] );
    // std::cout << tmp_log << std::endl;

    sum_log_ratio += tmp_log;
  }
  KL_est = (double)dimX/(double)xN * sum_log_ratio + log( (double)yN / ((double)xN-1.0 ) );

  // Deallocate memory
  annDeallocPts( dataX );
  annDeallocPts( dataY );
  delete [] distsX;
  delete [] distsXY;

  return KL_est;
}

#endif // QUESO_HAS_ANN

#endif // __UQ_INFO_THEORY_H__
