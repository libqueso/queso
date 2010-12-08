//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010 The PECOS Development Team
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
// $Id$
//
//--------------------------------------------------------------------------

#include <example_compute.h>
#include <uqGslMatrix.h>
#include <uqVectorRV.h>

void compute(const uqFullEnvironmentClass& env) {
  // Step 1 of 9: Instantiate the parameter space
  uqVectorSpaceClass<uqGslVectorClass,uqGslMatrixClass>
    paramSpace(env, "param_", 2, NULL);

  // Step 2 of 9: Instantiate the parameter domain
  uqGslVectorClass paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  uqGslVectorClass paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  uqBoxSubsetClass<uqGslVectorClass,uqGslMatrixClass>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  // Step 3 of 9: Instantiate the vector RV
  uqGslVectorClass meanVector(paramSpace.zeroVector());
  meanVector[0] = -1;
  meanVector[1] =  2;
  uqGslMatrixClass covMatrix = uqGslMatrixClass(paramSpace.zeroVector());
  covMatrix(0,0) = 4.; covMatrix(0,1) = 0.;
  covMatrix(1,0) = 0.; covMatrix(1,1) = 7.;
  uqGaussianVectorRVClass<uqGslVectorClass,uqGslMatrixClass>
    auxRv("", paramDomain,meanVector,covMatrix);

  // Step 4 of 9: Instantiate the vector sequence
  const unsigned int numSamples = 1000;

  uqSequenceOfVectorsClass<uqGslVectorClass,uqGslMatrixClass>
    auxSeq(paramSpace,numSamples,"aux_seq");

  // Step 5 of 9: Populate the vector sequence
  uqGslVectorClass auxVec(paramSpace.zeroVector());
  for (unsigned int i = 0; i < auxSeq.subSequenceSize(); ++i) {
    auxRv.realizer().realization(auxVec);
    auxSeq.setPositionValues(i,auxVec);
  }

#if 0
  // Step 6 of 9: Deal with BrooksGelman
  std::vector<double> convMeasure(numSamples,0.0);
  for (unsigned int i = 0; i < auxSeq.subSequenceSize(); ++i) {
    if ((i >= 1) && (env.numSubEnvironments() > 1)) {
      convMeasure[i] = auxSeq.estimateConvBrooksGelman(0,i);
    }
  }
  std::set<unsigned int> auxSet;
  auxSet.insert(0);
  auxSet.insert(1);
  auxSeq.subWriteContents("anyname",auxSet);
  if (env.inter0Rank() == 0) {
    std::ofstream dataout( "convergence.m", std::ios::out );
    dataout << "clear all" << std::endl;
    dataout << "close all" << std::endl;
    dataout << "index = [ " << 1 << std::endl;
    for (int i = 2; i < numSamples-1; i++) {
      dataout << i << std::endl;
    }
    dataout << numSamples-1 << " ];" << std::endl;
    dataout << " data = [ " << convMeasure[1] << std::endl;
    for (int i = 2; i < numSamples-1; i++) {
      dataout << convMeasure[i] << std::endl;
    }
    dataout << convMeasure[ numSamples-1] << "];" << std::endl;
    dataout << "plot( index, data, 'b-', 'LineWidth', 2 )" << std::endl;
    dataout << "xlabel( 'Iteration', 'FontSize', 16 )" << std::endl;
    dataout << "ylabel( 'BG-Convergence', 'FontSize', 16 )" << std::endl;
    dataout << "title( 'Brooks-Gelman Convergence, Gaussian RV, " 
	    << env.numSubEnvironments() << " Sequences', 'FontSize', 16 )" << std::endl;
    dataout << "print -depsc BGConv" << env.numSubEnvironments() << ".eps" << std::endl;
    std::cout <<"convMeasure = " << convMeasure[numSamples-1] << std::endl;
  }
#endif
  
  // Step 7 of 9: Compute min, max, mean, covariance and correlation matrices
  uqGslVectorClass minVec (paramSpace.zeroVector());
  uqGslVectorClass maxVec (paramSpace.zeroVector());
  auxSeq.unifiedMinMax(0,auxSeq.subSequenceSize(),minVec,maxVec);

  uqGslVectorClass meanVec(paramSpace.zeroVector());
  auxSeq.unifiedMean(0,auxSeq.subSequenceSize(),meanVec);

  uqGslMatrixClass covarianceMatrix  = uqGslMatrixClass(paramSpace.zeroVector());
  uqGslMatrixClass correlationMatrix = uqGslMatrixClass(paramSpace.zeroVector());
  uqComputeCovCorrMatricesBetweenVectorSequences(auxSeq,
                                                 auxSeq,
                                                 auxSeq.subSequenceSize(),
                                                 covarianceMatrix,
                                                 correlationMatrix);

  if (env.fullRank() == 0) {
    std::cout << "\n minVec = "  << minVec
              << "\n maxVec = "  << maxVec
              << "\n meanVec = " << meanVec
              << "\n covMat = "  << covarianceMatrix
              << "\n corrMat = " << correlationMatrix
              << std::endl;
  }
  
  // Step 8 of 9: Compute cdf accuracy
  std::vector<uqGslVectorClass*> cdfStaccVecs   (numSamples,NULL);
  std::vector<uqGslVectorClass*> cdfStaccVecsUp (numSamples,NULL);
  std::vector<uqGslVectorClass*> cdfStaccVecsLow(numSamples,NULL);
  std::vector<uqGslVectorClass*> sortedDataVecs (numSamples,NULL);
  std::cout << "Calling subCdfStacc()..."
            << std::endl;
  struct timeval timevalBegin;
  gettimeofday(&timevalBegin, NULL);
  auxSeq.subCdfStacc(0,
                     cdfStaccVecs,
                     cdfStaccVecsUp,
                     cdfStaccVecsLow,
                     sortedDataVecs);
  double traceT = uqMiscGetEllapsedSeconds(&timevalBegin);
  std::cout << "Returned from subCdfStacc()"
            << " after " << traceT << " seconds"
            << std::endl;

  if (env.fullRank() == 0) {
    std::ofstream outdata1;
    std::ofstream outdata2;
    std::ofstream outdata3;
    std::ofstream outdata4;
    outdata1.open("./cdf.dat"       );
    outdata2.open("./cdfup.dat"     );
    outdata3.open("./cdfdown.dat"   );
    outdata4.open("./sortedData.dat");
    for (unsigned int i = 0; i < numSamples; ++i) {
      //std::cout << *cdfStaccVecs   [i] << std::endl;
      //std::cout << *cdfStaccVecsUp [i] << std::endl;
      //std::cout << *cdfStaccVecsLow[i] << std::endl;
      //std::cout << *sortedDataVecs [i] << std::endl;
      outdata1 << *cdfStaccVecs   [i] << std::endl;
      outdata2 << *cdfStaccVecsUp [i] << std::endl;
      outdata3 << *cdfStaccVecsLow[i] << std::endl;
      outdata4 << *sortedDataVecs [i] << std::endl;
    }
    outdata1.close();
    outdata2.close();
    outdata3.close();
    outdata4.close();
    //std::cout << "\n cdfmean = " << cdfStaccVecs << std::endl;
  }

  for (unsigned int i = 0; i < numSamples; ++i) {
    delete cdfStaccVecs   [i];
    delete cdfStaccVecsUp [i];
    delete cdfStaccVecsLow[i];
    delete sortedDataVecs [i];
  }

#if 0
  uqGslVectorClass deltaVec(maxVec-minVec);
  deltaVec *= (1./(double) (auxSize-1));
  std::vector<uqGslVectorClass*> evalPositionsVecs(auxSize,NULL);
  for (unsigned int i = 0; i < auxSize; ++i) {
    evalPositionsVecs[i] = new uqGslVectorClass(paramSpace.zeroVector());
    *(evalPositionsVecs[i]) = minVec + ((double) i)*deltaVec;
  }
  for (unsigned int i = 0; i < auxSize; ++i) {
    delete evalPositionsVecs[i];
  }
#endif

  //Return
  return;
}
