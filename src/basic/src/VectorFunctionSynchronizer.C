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

#include <queso/VectorFunctionSynchronizer.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor -----------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>::VectorFunctionSynchronizer(
  const BaseVectorFunction<P_V,P_M,Q_V,Q_M>& inputFunction,
  const P_V&                                        auxPVec,
  const Q_V&                                        auxQVec)
  :
  m_env           (inputFunction.domainSet().env()),
  m_vectorFunction(inputFunction),
  m_auxPVec       (auxPVec),
  m_auxQVec       (auxQVec)
{
}

// Destructor ---------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>::~VectorFunctionSynchronizer()
{
}
// Math methods -------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
const VectorSet<P_V,P_M>&
VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>::domainSet() const
{
  return m_vectorFunction.domainSet();
}

// Sync methods -------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
void
VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>::callFunction(
  const P_V*                    vecValues,
  const P_V*                    vecDirection,
        Q_V*                    imageVector,
        DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
        DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
        DistArray<P_V*>* hessianEffects) const
{
  if ((m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) &&
      (m_auxPVec.numOfProcsForStorage() == 1                                 ) &&
      (m_auxQVec.numOfProcsForStorage() == 1                                 )) {
    bool stayInRoutine = true;
    do {
      const P_V*                    internalValues    = NULL;
      const P_V*                    internalDirection = NULL;
            Q_V*                    internalImageVec  = NULL;
            DistArray<P_V*>* internalGrads     = NULL; // Yes, 'P_V'
            DistArray<P_M*>* internalHessians  = NULL; // Yes, 'P_M'
            DistArray<P_V*>* internalEffects   = NULL;

      /////////////////////////////////////////////////
      // Broadcast 1 of 3
      /////////////////////////////////////////////////
      // bufferChar[0] = '0' or '1' (vecValues       is NULL or not)
      // bufferChar[1] = '0' or '1' (vecDirection    is NULL or not)
      // bufferChar[2] = '0' or '1' (imageVector     is NULL or not)
      // bufferChar[3] = '0' or '1' (gradVectors     is NULL or not)
      // bufferChar[4] = '0' or '1' (hessianMatrices is NULL or not)
      // bufferChar[5] = '0' or '1' (hessianEffects  is NULL or not)
      std::vector<char> bufferChar(6,'0');

      if (m_env.subRank() == 0) {
        if ((vecValues != NULL)) queso_require_msg(imageVector, "imageVector should not be NULL");
        internalValues    = vecValues;
        internalDirection = vecDirection;
        internalImageVec  = imageVector;
        internalGrads     = gradVectors;
        internalHessians  = hessianMatrices;
        internalEffects   = hessianEffects;

        if (internalValues    != NULL) bufferChar[0] = '1';
        if (internalDirection != NULL) bufferChar[1] = '1';
        if (internalImageVec  != NULL) bufferChar[2] = '1';
        if (internalGrads     != NULL) bufferChar[3] = '1';
        if (internalHessians  != NULL) bufferChar[4] = '1';
        if (internalEffects   != NULL) bufferChar[5] = '1';
      }

      m_env.subComm().syncPrintDebugMsg("In VectorFunctionSynchronizer<V,M>::callFunction(), just before char Bcast()",3,3000000);
      //if (m_env.subId() != 0) while (true) sleep(1);

      int count = (int) bufferChar.size();
      m_env.subComm().Bcast((void *) &bufferChar[0], count, RawValue_MPI_CHAR, 0,
                            "VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>::callFunction()",
                            "failed broadcast 1 of 3");

      if (bufferChar[0] == '1') {
        ///////////////////////////////////////////////
        // Broadcast 2 of 3
        ///////////////////////////////////////////////

        // bufferDouble[0...] = contents for (eventual) vecValues
        std::vector<double> bufferDouble(m_auxPVec.sizeLocal(),0);

        if (m_env.subRank() == 0) {
          for (unsigned int i = 0; i < internalValues->sizeLocal(); ++i) {
            bufferDouble[i] = (*internalValues)[i];
          }
        }

        count = (int) bufferDouble.size();
        m_env.subComm().Bcast((void *) &bufferDouble[0], count, RawValue_MPI_DOUBLE, 0,
                              "VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>::callFunction()",
                              "failed broadcast 2 of 3");

        if (m_env.subRank() != 0) {
          P_V tmpPVec(m_auxPVec);
          for (unsigned int i = 0; i < tmpPVec.sizeLocal(); ++i) {
            tmpPVec[i] = bufferDouble[i];
          }
          internalValues = new P_V(tmpPVec);
        }

        if (bufferChar[1] == '1') {
          /////////////////////////////////////////////
          // Broadcast 3 of 3
          /////////////////////////////////////////////
          // bufferDouble[0...] = contents for (eventual) vecDirection

          if (m_env.subRank() == 0) {
            for (unsigned int i = 0; i < internalDirection->sizeLocal(); ++i) {
              bufferDouble[i] = (*internalDirection)[i];
            }
          }

          count = (int) bufferDouble.size();
          m_env.subComm().Bcast((void *) &bufferDouble[0], count, RawValue_MPI_DOUBLE, 0,
                                "VectorFunctionSynchronizer<P_V,P_M,Q_V,Q_M>::callFunction()",
                                "failed broadcast 3 of 3");

          if (m_env.subRank() != 0) {
            P_V tmpPVec(m_auxPVec);
            for (unsigned int i = 0; i < tmpPVec.sizeLocal(); ++i) {
              tmpPVec[i] = bufferDouble[i];
            }
            internalDirection = new P_V(tmpPVec);
          }
        }

        ///////////////////////////////////////////////
        // All processors now call 'vectorFunction()'
        ///////////////////////////////////////////////
        if (m_env.subRank() != 0) {
          if (bufferChar[2] == '1') internalImageVec = new Q_V(m_auxQVec);
        //if (bufferChar[3] == '1') internalGrads    = new P_V(m_auxPVec);
        //if (bufferChar[4] == '1') internalHessians = new P_M(m_auxPVec);
        //if (bufferChar[5] == '1') internalEffects  = new P_V(m_auxPVec);
        }

        m_env.subComm().Barrier();
        m_vectorFunction.compute(*internalValues,
                                 internalDirection,
                                 *internalImageVec,
                                 internalGrads,
                                 internalHessians,
                                 internalEffects);
      }

      /////////////////////////////////////////////////
      // Prepare to exit routine or to stay in it
      /////////////////////////////////////////////////
      if (m_env.subRank() == 0) {
        stayInRoutine = false; // Always for processor 0
      }
      else {
        if (internalValues    != NULL) delete internalValues;
        if (internalDirection != NULL) delete internalDirection;
        if (internalImageVec  != NULL) delete internalImageVec;
      //if (internalGrads     != NULL) delete internalGrads;
      //if (internalHessians  != NULL) delete internalHessians;
      //if (internalEffects   != NULL) delete internalEffects;

        stayInRoutine = (vecValues == NULL) && (bufferChar[0] == '1');
      }
    } while (stayInRoutine);
  }
  else {
    queso_require_msg(!((vecValues == NULL) || (imageVector == NULL)), "Neither vecValues nor imageVector should not be NULL");
    queso_require_equal_to_msg(m_auxPVec.numOfProcsForStorage(), m_auxQVec.numOfProcsForStorage(), "Number of processors required for storage should be the same");

    m_env.subComm().Barrier();
    m_vectorFunction.compute(*vecValues,
                             vecDirection,
                             *imageVector,
                             gradVectors,
                             hessianMatrices,
                             hessianEffects);
  }

  return;
}

}  // End namespace QUESO

template class QUESO::VectorFunctionSynchronizer<QUESO::GslVector, QUESO::GslMatrix, QUESO::GslVector, QUESO::GslMatrix>;
