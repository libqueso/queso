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

#include <queso/VectorRV.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

#include <queso/Defines.h>
#include <gsl/gsl_sf_psi.h> // todo: take specificity of gsl_, i.e., make it general (gsl or boost or etc)
#include <queso/InfoTheory_impl.h>

namespace QUESO {

// Default constructor -----------------------------
template<class V, class M>
BaseVectorRV<V,M>::BaseVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>& imageSet)
  :
  m_env       (imageSet.env()),
  m_prefix    ((std::string)(prefix)+"rv_"),
  m_imageSet  (imageSet),
  m_pdf       (NULL),
  m_realizer  (NULL),
  m_subCdf    (NULL),
  m_unifiedCdf(NULL),
  m_mdf       (NULL)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering BaseVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving BaseVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
BaseVectorRV<V,M>::~BaseVectorRV()
{
  //if (m_mdf       ) delete m_mdf;
  //if (m_subCdf    ) delete m_subCdf;
  //if (m_unifiedCdf) delete m_unifiedCdf;
  //if (m_realizer  ) delete m_realizer;
  //if (m_pdf       ) delete m_pdf;
}
// RV handling-methods ------------------------------

template <class V, class M>
const BaseEnvironment&
BaseVectorRV<V,M>::env() const
{
  return m_env;
}
//---------------------------------------------------
template<class V, class M>
const VectorSet<V,M>&
BaseVectorRV<V,M>::imageSet() const
{
  return m_imageSet;
}
//---------------------------------------------------
template<class V, class M>
const BaseJointPdf<V,M>&
BaseVectorRV<V,M>::pdf() const
{
  queso_require_msg(m_pdf, "m_pdf is NULL");

  return *m_pdf;
}
//---------------------------------------------------
template<class V, class M>
bool
BaseVectorRV<V,M>::has_realizer() const
{
  return (m_realizer != NULL);
}
//---------------------------------------------------
template<class V, class M>
const BaseVectorRealizer<V,M>&
BaseVectorRV<V,M>::realizer() const
{
  queso_require_msg(m_realizer, "m_realizer is NULL");

  return *m_realizer;
}
//---------------------------------------------------
template<class V, class M>
const BaseVectorCdf<V,M>&
BaseVectorRV<V,M>::subCdf() const
{
  queso_require_msg(m_subCdf, "m_subCdf is NULL");

  return *m_subCdf;
}
//---------------------------------------------------
template<class V, class M>
const BaseVectorCdf<V,M>&
BaseVectorRV<V,M>::unifiedCdf() const
{
  queso_require_msg(m_unifiedCdf, "m_unifiedCdf is NULL");

  return *m_unifiedCdf;
}
//---------------------------------------------------
template<class V, class M>
const BaseVectorMdf<V,M>&
BaseVectorRV<V,M>::mdf() const
{
  queso_require_msg(m_mdf, "m_mdf is NULL");

  return *m_mdf;
}

//---------------------------------------------------
#ifdef QUESO_HAS_ANN
template <class V, class M>
double
BaseVectorRV<V,M>::estimateENT_ANN() const
{
  ANNpointArray data;
  double* dists;
  double ENT_est;

  // FIXME: these default values should be stored in the
  // QUESO input file ( create a InfoTheoryOptions )
  unsigned int k = UQ_INFTH_ANN_KNN;
  double eps = UQ_INFTH_ANN_EPS;

  // here it is assumed that the entropy for the
  // entire joint RV will be computed
  unsigned int dim = this->imageSet().vectorSpace().dimGlobal();

  // FIXME: get the number already stored, otherwise
  // use the default value
  unsigned int N = this->realizer().subPeriod();
  if( N == 0 ) {
    N = UQ_INFTH_ANN_NO_SMP;
  }

  // allocate memory
  data = annAllocPts(N,dim);
  dists = new double[N];

  // copy samples in the ANN data structure
  V smpRV( this->imageSet().vectorSpace().zeroVector() );
  for( unsigned int i = 0; i < N; i++ ) {
    // get a sample from the distribution
    this->realizer().realization( smpRV );
    // copy the vector values in the ANN data structure
    for( unsigned int j = 0; j < dim; j++ ) {
      data[ i ][ j ] = smpRV[ j ];
    }
  }

  // get distance to knn for each point
  // (k+1) because the 1st nn is itself
  distANN_XY( data, data, dists, dim, dim, N, N, k+1, eps );

  // compute the entropy estimate using the L-infinity (Max) norm
  // so no need for the adjustment of the mass of the hyperball
  // this has to be enforced before compiling the ANN lib
  double sum_log_dist = 0.0;
  for( unsigned int i = 0; i < N; i++ ) {
    if( dists[ i ] > 0 ) {
      sum_log_dist += log( 2.0*dists[ i ] );
    }
  }
  ENT_est = - gsl_sf_psi_int( k ) + gsl_sf_psi_int( N ) + (double)dim / (double)N * sum_log_dist; // todo: take specificity of gsl_, i.e., make it general (gsl or boost or etc)

  // deallocate memory
  delete [] dists;
  annDeallocPts( data );

  return ENT_est;
}
#endif // QUESO_HAS_ANN

//---------------------------------------------------
// Method declared outside class definition ---------
//---------------------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
void
ComputeCovCorrMatricesBetweenVectorRvs(
  const BaseVectorRV<P_V,P_M>& paramRv,
  const BaseVectorRV<Q_V,Q_M>& qoiRv,
        unsigned int                  localNumSamples,
        P_M&                          pqCovMatrix,
        P_M&                          pqCorrMatrix)
{
  // Check input data consistency
  bool useOnlyInter0Comm = (paramRv.imageSet().vectorSpace().numOfProcsForStorage() == 1) &&
                           (qoiRv.imageSet().vectorSpace().numOfProcsForStorage()   == 1);

  queso_require_msg(useOnlyInter0Comm, "parallel vectors not supported yet");

  unsigned int numRows = paramRv.imageSet().vectorSpace().dim();
  unsigned int numCols = qoiRv.imageSet().vectorSpace().dim();

  queso_require_msg(!((numRows != pqCovMatrix.numRows()) || (numCols != pqCovMatrix.numCols())), "inconsistent dimensions for covariance matrix");

  queso_require_msg(!((numRows != pqCorrMatrix.numRows()) || (numCols != pqCorrMatrix.numCols())), "inconsistent dimensions for correlation matrix");

  queso_require_msg(!((localNumSamples > paramRv.realizer().period()) || (localNumSamples > qoiRv.realizer().period())), "localNumSamples is too large");

  // For both P and Q vector sequences: fill them
  P_V tmpP(paramRv.imageSet().vectorSpace().zeroVector());
  Q_V tmpQ(qoiRv.imageSet().vectorSpace().zeroVector());

  SequenceOfVectors<P_V,P_M> localWorkingPSeq(paramRv.imageSet().vectorSpace(),
                                                     localNumSamples,
                                                     "covTmpP");
  SequenceOfVectors<Q_V,Q_M> localWorkingQSeq(qoiRv.imageSet().vectorSpace(),
                                                     localNumSamples,
                                                     "covTmpQ");
  for (unsigned int k = 0; k < localNumSamples; ++k) {
    paramRv.realizer().realization(tmpP);
    localWorkingPSeq.setPositionValues(k,tmpP);

    qoiRv.realizer().realization(tmpQ);
    localWorkingQSeq.setPositionValues(k,tmpQ);
  }

  ComputeCovCorrMatricesBetweenVectorSequences(localWorkingPSeq,
                                                 localWorkingQSeq,
                                                 localNumSamples,
                                                 pqCovMatrix,
                                                 pqCorrMatrix);

  return;
}

}  // End namespace QUESO

template class QUESO::BaseVectorRV<QUESO::GslVector,QUESO::GslMatrix>;
