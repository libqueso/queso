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
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_VECTOR_FUNCTION_SYNCHRONIZER_H__
#define __UQ_VECTOR_FUNCTION_SYNCHRONIZER_H__

namespace QUESO {

/*! \file uqVectorFunctionSynchronizer.h
 * \brief Class for synchronizing the calls of vector-valued functions
 * 
 * \class uqVectorFunctionSynchronizerClass
 * \brief A templated class for synchronizing the calls of vector-valued functions.
 *
 * This class creates a synchronization point among processes which call vector-valued 
 * functions. This means that all processes must reach a point in their code before they 
 * can all begin executing again. */

template <class P_V, class P_M, class Q_V, class Q_M>
class uqVectorFunctionSynchronizerClass
{
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  uqVectorFunctionSynchronizerClass(const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& inputFunction,
                                    const P_V&                                        auxPVec,
                                    const Q_V&                                        auxQVec);
  //! Destructor
 ~uqVectorFunctionSynchronizerClass();
  //@}
 
  //! @name Mathematical methods
  //@{  
  //! Access to the domain set of the vector-valued function which will be synchronized.
  const uqVectorSetClass<P_V,P_M>& domainSet() const;
  //@}
  
  //! @name Sync method
  //@{  
  //! Calls the vector-valued function which will be synchronized.
  /*! This procedure  forms a barrier, and no processes in the communicator can pass the 
   * barrier until all of them call the function. */
  void callFunction(const P_V*                    vecValues,
                    const P_V*                    vecDirection,
                          Q_V*                    imageVector,
                          uqDistArrayClass<P_V*>* gradVectors,     // Yes, 'P_V'
                          uqDistArrayClass<P_M*>* hessianMatrices, // Yes, 'P_M'
                          uqDistArrayClass<P_V*>* hessianEffects) const;
  //@}
private:
  const uqBaseEnvironmentClass&                     m_env;
  const uqBaseVectorFunctionClass<P_V,P_M,Q_V,Q_M>& m_vectorFunction;
  const P_V&                                        m_auxPVec;
  const Q_V&                                        m_auxQVec;
};
// Default constructor -----------------------------
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
// Destructor ---------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::~uqVectorFunctionSynchronizerClass()
{
}
// Math methods -------------------------------------
template<class P_V, class P_M, class Q_V, class Q_M>
const uqVectorSetClass<P_V,P_M>&
uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::domainSet() const
{
  return m_vectorFunction.domainSet();
}
// Sync methods -------------------------------------
template <class P_V, class P_M, class Q_V, class Q_M>
void
uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction(
  const P_V*                    vecValues,
  const P_V*                    vecDirection,
        Q_V*                    imageVector,
        uqDistArrayClass<P_V*>* gradVectors,     // Yes, 'P_V'
        uqDistArrayClass<P_M*>* hessianMatrices, // Yes, 'P_M'
        uqDistArrayClass<P_V*>* hessianEffects) const
{
  if ((m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) &&
      (m_auxPVec.numOfProcsForStorage() == 1                                 ) &&
      (m_auxQVec.numOfProcsForStorage() == 1                                 )) {
    bool stayInRoutine = true;
    do {
      const P_V*                    internalValues    = NULL;
      const P_V*                    internalDirection = NULL;
            Q_V*                    internalImageVec  = NULL;
            uqDistArrayClass<P_V*>* internalGrads     = NULL; // Yes, 'P_V'
            uqDistArrayClass<P_M*>* internalHessians  = NULL; // Yes, 'P_M'
            uqDistArrayClass<P_V*>* internalEffects   = NULL;

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
                            m_env.worldRank(),
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

      m_env.subComm().syncPrintDebugMsg("In uqVectorFunctionSynchronizerClass<V,M>::callFunction(), just before char Bcast()",3,3000000);
      //if (m_env.subId() != 0) while (true) sleep(1);

      int count = (int) bufferChar.size();
      m_env.subComm().Bcast((void *) &bufferChar[0], count, uqRawValue_MPI_CHAR, 0,
                            "uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction()",
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
        m_env.subComm().Bcast((void *) &bufferDouble[0], count, uqRawValue_MPI_DOUBLE, 0,
                              "uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction()",
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
          m_env.subComm().Bcast((void *) &bufferDouble[0], count, uqRawValue_MPI_DOUBLE, 0,
                                "uqVectorFunctionSynchronizerClass<P_V,P_M,Q_V,Q_M>::callFunction()",
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
    UQ_FATAL_TEST_MACRO((vecValues == NULL) || (imageVector == NULL),
                        m_env.worldRank(),
                        "uqVectorFunctionSynchronizerClass<V,M>::callFunction()",
                        "Neither vecValues nor imageVector should not be NULL");
    UQ_FATAL_TEST_MACRO((m_auxPVec.numOfProcsForStorage() != m_auxQVec.numOfProcsForStorage()),
                        m_env.worldRank(),
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

}  // End namespace QUESO

#endif // __UQ_VECTOR_FUNCTION_SYNCHRONIZER_H__
