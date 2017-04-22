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

#include <queso/Defines.h>

#ifdef QUESO_HAS_ANN
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GaussianVectorRV.h>
#include <queso/InfoTheory.h>

#include <gsl/gsl_sf_psi.h> // todo: take specificity of gsl_, i.e., make it general (gsl or boost or etc)

namespace QUESO {

//*****************************************************
// Function: distANN_XY
//*****************************************************
void distANN_XY( const ANNpointArray dataX, const ANNpointArray dataY,
		 double* distsXY,
		 unsigned int dimX, unsigned int dimY,
		 unsigned int xN, unsigned int yN,
		 unsigned int k, double eps )
{

  ANNkd_tree* kdTree;
  ANNdistArray nnDist;
  ANNidxArray nnIdx;

  // Allocate memory
  nnIdx = new ANNidx[ k+1 ];
  nnDist = new ANNdist[ k+1 ];
  kdTree = new ANNkd_tree( dataY, yN, dimY );

  // Get the distances to all the points
  for( unsigned int i = 0; i < xN ; i++ )
    {
      kdTree->annkSearch( dataX[ i ], k+1, nnIdx, nnDist, eps );

      double my_dist = nnDist[ k ];

      // check to see if the dist is zero (query point same as the kNN)
      // if so find the next k that gives the next positive distance
      if( my_dist == 0.0 )
	{
	  ANNdistArray nnDist_tmp = new ANNdist[ yN ];
	  ANNidxArray nnIdx_tmp = new ANNidx[ yN ];
	  kdTree->annkSearch( dataX[ i ], yN, nnIdx_tmp, nnDist_tmp, eps );

	  for( unsigned int my_k = k + 1; my_k < yN; ++my_k )
	    if( nnDist_tmp[ my_k ] > 0.0 )
	      {
		my_dist = nnDist_tmp[ my_k ];
		break;
	      }
	  delete [] nnIdx_tmp;
	  delete [] nnDist_tmp;
	}

      distsXY[ i ] = my_dist;
    }

  // Deallocate memory
  delete [] nnIdx;
  delete [] nnDist;
  delete kdTree;
  annClose();

  return;
}

//*****************************************************
// Function: normalizeANN_XY
// (used by Mutual Information - marginal normalization
//  it does not destroy the correlations)
//*****************************************************
void normalizeANN_XY( ANNpointArray dataXY, unsigned int dimXY,
		      ANNpointArray dataX, unsigned int dimX,
		      ANNpointArray dataY, unsigned int dimY,
		      unsigned int N )
{

  ANNpoint meanXY, stdXY;

  meanXY = annAllocPt(dimXY); // is init with 0
  stdXY = annAllocPt(dimXY); // is init with 0

  // get mean
  for( unsigned int i = 0; i < N; i++ ) {
    for( unsigned int j = 0; j < dimX; j++ ) {
      meanXY[ j ] += dataXY[ i ][ j ];
    }
    for( unsigned int j = 0; j < dimY; j++ ) {
      meanXY[ dimX + j ] += dataXY[ i ][ dimX + j ];
    }
  }
  for( unsigned int j = 0; j < dimXY; j++ ) {
    meanXY[ j ] = meanXY[ j ] / (double)N;
  }

  // get standard deviation
  for( unsigned int i = 0; i < N; i++ ) {
    for( unsigned int j = 0; j < dimXY; j++ ) {
      stdXY[ j ] += pow( dataXY[ i ][ j ] - meanXY[ j ], 2.0 );
    }
  }
  for( unsigned int j = 0; j < dimXY; j++ ) {
    stdXY[ j ] = sqrt( stdXY[ j ] / ((double)N-1.0) );
  }

  /*
  std::cout << "Mean RV - ";
  annPrintPt( meanXY, dimXY, std::cout );
  std::cout << std::endl << "Std RV - ";
  annPrintPt( stdXY, dimXY, std::cout );
  std::cout << std::endl;
  */

  // get normalized points and populate marginals
  for( unsigned int i = 0; i < N; i++ ) {
    // normalize joint
    for( unsigned int j = 0; j < dimXY; j++ ) {
      dataXY[ i ][ j ] = ( dataXY[ i ][ j ] - meanXY[ j ] ) / stdXY[ j ];
    }

    // populate marginals
    for( unsigned int j = 0; j < dimX; j++ ) {
      dataX[ i ][ j ] = dataXY[ i ][ j ];
    }
    for( unsigned int j = 0; j < dimY; j++ ) {
      dataY[ i ][ j ] = dataXY[ i ][ dimX + j ];
    }
  }

}

//*****************************************************
// Function: computeMI_ANN
//*****************************************************
double computeMI_ANN( ANNpointArray dataXY,
		      unsigned int dimX, unsigned int dimY,
		      unsigned int k, unsigned int N, double eps )
{

  ANNpointArray dataX, dataY;
  double* distsXY;
  double MI_est;
  ANNkd_tree* kdTreeX;
  ANNkd_tree* kdTreeY;

  unsigned int dimXY = dimX + dimY;

  // Allocate memory
  dataX = annAllocPts(N,dimX);
  dataY = annAllocPts(N,dimY);
  distsXY = new double[N];

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

  return MI_est;

}

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
    queso_error_msg("Error-KL: the dimensions should agree");
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

  // sanity check
  if( dimX != dimY ) {
    queso_error_msg("Error-CE: the dimensions should agree");
  }

  // Allocate memory
  dataX = annAllocPts( xN, dimX );
  dataY = annAllocPts( yN, dimY );
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
  distANN_XY( dataX, dataY, distsXY, dimX, dimY, xN, yN, k, eps );

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
