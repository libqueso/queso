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

#ifndef __UQ_VECTOR_FUNCTION_SYNCHRONIZER_H__
#define __UQ_VECTOR_FUNCTION_SYNCHRONIZER_H__

template <class P_V, class P_M, class Q_V, class Q_M>
class uqVectorFunctionSynchronizerClass
{
public:
  uqVectorFunctionSynchronizerClass(const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& inputFunction,
                                    const P_V&                                        auxPVec,
                                    const Q_V&                                        auxQVec);
 ~uqVectorFunctionSynchronizerClass();

  const uqVectorSetClass<P_V,P_M>& domainSet() const;
  void callFunction(const P_V*                        vecValues,
                    const P_V*                        vecDirection,
                          Q_V*                        imageVector,
                          EpetraExt::DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
                          EpetraExt::DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
                          EpetraExt::DistArray<P_V*>* hessianEffects) const;
private:
  const uqBaseEnvironmentClass&                     m_env;
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& m_vectorFunction;
  const P_V&                                        m_auxPVec;
  const Q_V&                                        m_auxQVec;
};

template <class P_V, class P_M, class Q_V, class Q_M>
uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::uqVectorFunctionSynchronizerClass(
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& inputFunction,
  const P_V&                                        auxPVec,
  const Q_V&                                        auxQVec)
  :
  m_env           (inputFunction.domainSet().env()),
  m_vectorFunction(inputFunction),
  m_auxPVec       (auxPVec),
  m_auxQVec       (auxQVec)
{
}

template <class P_V, class P_M, class Q_V, class Q_M>
uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::~uqVectorFunctionSynchronizerClass()
{
}

template<class P_V, class P_M, class Q_V, class Q_M>
const uqVectorSetClass<P_V,P_M>&
uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::domainSet() const
{
  return m_vectorFunction.domainSet();
}

template <class P_V, class P_M, class Q_V, class Q_M>
void
uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction(
  const P_V*                        vecValues,
  const P_V*                        vecDirection,
        Q_V*                        imageVector,
        EpetraExt::DistArray<P_V*>* gradVectors,     // Yes, 'P_V'
        EpetraExt::DistArray<P_M*>* hessianMatrices, // Yes, 'P_M'
        EpetraExt::DistArray<P_V*>* hessianEffects) const
{
  if ((m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) &&
      (m_auxPVec.numberOfProcessorsRequiredForStorage() == 1                 ) &&
      (m_auxQVec.numberOfProcessorsRequiredForStorage() == 1                 )) {
    bool stayInRoutine = true;
    do {
      const P_V*                        internalValues    = NULL;
      const P_V*                        internalDirection = NULL;
            Q_V*                        internalImageVec  = NULL;
            EpetraExt::DistArray<P_V*>* internalGrads     = NULL; // Yes, 'P_V'
            EpetraExt::DistArray<P_M*>* internalHessians  = NULL; // Yes, 'P_M'
            EpetraExt::DistArray<P_V*>* internalEffects   = NULL;

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
        UQ_FATAL_TEST_MACRO((vecValues != NULL) && (imageVector == NULL),
                            m_env.rank(),
                            "uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction()",
                            "imageVector should not be NULL");
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

      m_env.syncPrintDebugMsg("In uqVectorFunctionSynchronizerClass<V,M>::callFunction(), just before char Bcast()",3,3000000,m_env.subComm());
      //if (m_env.subId() != 0) while (true) sleep(1);

      int count = (int) bufferChar.size();
      int mpiRC = MPI_Bcast ((void *) &bufferChar[0], count, MPI_CHAR, 0, m_env.subComm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.rank(),
                          "uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction()",
                          "failed broadcast 1 of 3");

      if (bufferChar[0] == '1') {
        ///////////////////////////////////////////////
        // Broadcast 2 of 3
        ///////////////////////////////////////////////

        // bufferDouble[0...] = contents for (eventual) vecValues
        std::vector<double> bufferDouble(m_auxPVec.size(),0);

        if (m_env.subRank() == 0) {
          for (unsigned int i = 0; i < internalValues->size(); ++i) {
            bufferDouble[i] = (*internalValues)[i];
          }
        }

        count = (int) bufferDouble.size();
        mpiRC = MPI_Bcast ((void *) &bufferDouble[0], count, MPI_DOUBLE, 0, m_env.subComm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.rank(),
                            "uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction()",
                            "failed broadcast 2 of 3");

        if (m_env.subRank() != 0) {
          P_V tmpPVec(m_auxPVec);
          for (unsigned int i = 0; i < tmpPVec.size(); ++i) {
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
            for (unsigned int i = 0; i < internalDirection->size(); ++i) {
              bufferDouble[i] = (*internalDirection)[i];
            }
          }

          count = (int) bufferDouble.size();
          mpiRC = MPI_Bcast ((void *) &bufferDouble[0], count, MPI_DOUBLE, 0, m_env.subComm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.rank(),
                              "uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction()",
                              "failed broadcast 3 of 3");

          if (m_env.subRank() != 0) {
            P_V tmpPVec(m_auxPVec);
            for (unsigned int i = 0; i < tmpPVec.size(); ++i) {
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
    UQ_FATAL_TEST_MACRO((vecValues == NULL) || (imageVector == NULL),
                        m_env.rank(),
                        "uqVectorFunctionSynchronizerClass<V,M>::callFunction()",
                        "Neither vecValues nor imageVector should not be NULL");
    UQ_FATAL_TEST_MACRO((m_auxPVec.numberOfProcessorsRequiredForStorage() != m_auxQVec.numberOfProcessorsRequiredForStorage()),
                        m_env.rank(),
                        "uqVectorFunctionSynchronizerClass<V,M>::callFunction()",
                        "Number of processors required for storage should be the same");

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

#endif // __UQ_VECTOR_FUNCTION_SYNCHRONIZER_H__
