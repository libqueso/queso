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

#ifndef __UQ_SCALAR_FUNCTION_SYNCHRONIZER_H__
#define __UQ_SCALAR_FUNCTION_SYNCHRONIZER_H__

template <class V, class M>
class uqScalarFunctionSynchronizerClass
{
public:
  uqScalarFunctionSynchronizerClass(const uqBaseScalarFunctionClass<V,M>& inputFunction,
                                    const V&                              auxVec);
 ~uqScalarFunctionSynchronizerClass();

  const uqVectorSetClass<V,M>& domainSet() const;
  double callFunction(const V* vecValues,
                      const V* vecDirection,
                            V* gradVector,
                            M* hessianMatrix,
                            V* hessianEffect) const;
private:
  const uqBaseEnvironmentClass&         m_env;
  const uqBaseScalarFunctionClass<V,M>& m_scalarFunction;
  const V&                              m_auxVec;
};

template <class V, class M>
uqScalarFunctionSynchronizerClass<V,M>::uqScalarFunctionSynchronizerClass(
  const uqBaseScalarFunctionClass<V,M>& inputFunction,
  const V&                              auxVec)
  :
  m_env           (inputFunction.domainSet().env()),
  m_scalarFunction(inputFunction),
  m_auxVec        (auxVec)
{
}

template <class V, class M>
uqScalarFunctionSynchronizerClass<V,M>::~uqScalarFunctionSynchronizerClass()
{
}

template<class V,class M>
const uqVectorSetClass<V,M>&
uqScalarFunctionSynchronizerClass<V,M>::domainSet() const
{
  return m_scalarFunction.domainSet();
}

template <class V,class M>
double
uqScalarFunctionSynchronizerClass<V,M>::callFunction(
  const V* vecValues,
  const V* vecDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  double result = 0.;

  if ((m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) &&
      (m_auxVec.numberOfProcessorsRequiredForStorage() == 1                  )) {
    bool stayInRoutine = true;
    do {
      const V* internalValues    = NULL;
      const V* internalDirection = NULL;
            V* internalGrad      = NULL;
            M* internalHessian   = NULL;
            V* internalEffect    = NULL;

      /////////////////////////////////////////////////
      // Broadcast 1 of 3
      /////////////////////////////////////////////////
      // bufferChar[0] = '0' or '1' (vecValues     is NULL or not)
      // bufferChar[1] = '0' or '1' (vecDirection  is NULL or not)
      // bufferChar[2] = '0' or '1' (gradVector    is NULL or not)
      // bufferChar[3] = '0' or '1' (hessianMatrix is NULL or not)
      // bufferChar[4] = '0' or '1' (hessianEffect is NULL or not)
      std::vector<char> bufferChar(5,'0');

      if (m_env.subRank() == 0) {
        internalValues    = vecValues;
        internalDirection = vecDirection;
        internalGrad      = gradVector;
        internalHessian   = hessianMatrix;
        internalEffect    = hessianEffect;

        if (internalValues    != NULL) bufferChar[0] = '1';
        if (internalDirection != NULL) bufferChar[1] = '1';
        if (internalGrad      != NULL) bufferChar[2] = '1';
        if (internalHessian   != NULL) bufferChar[3] = '1';
        if (internalEffect    != NULL) bufferChar[4] = '1';
      }

      m_env.syncPrintDebugMsg("In uqScalarFunctionSynchronizerClass<V,M>::callFunction(), just before char Bcast()",3,3000000,m_env.subComm());
      //if (m_env.subId() != 0) while (true) sleep(1);

      int count = (int) bufferChar.size();
      int mpiRC = MPI_Bcast ((void *) &bufferChar[0], count, MPI_CHAR, 0, m_env.subComm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          m_env.fullRank(),
                          "uqScalarFunctionSynchronizerClass<V,M>::callFunction()",
                          "failed broadcast 1 of 3");

      m_env.syncPrintDebugMsg("In uqScalarFunctionSynchronizerClass<V,M>::callFunction(), just after char Bcast()",3,3000000,m_env.subComm());
      //std::cout << "char contents = " << bufferChar[0] << " " << bufferChar[1] << " " << bufferChar[2] << " " << bufferChar[3] << " " << bufferChar[4]
      //          << std::endl;

      if (bufferChar[0] == '1') {
        ///////////////////////////////////////////////
        // Broadcast 2 of 3
        ///////////////////////////////////////////////

        // bufferDouble[0...] = contents for (eventual) vecValues
        std::vector<double> bufferDouble(m_auxVec.sizeLocal(),0.);

        if (m_env.subRank() == 0) {
          for (unsigned int i = 0; i < internalValues->sizeLocal(); ++i) {
            bufferDouble[i] = (*internalValues)[i];
          }
        }

        //m_env.fullComm().Barrier();
        //for (int i = 0; i < m_env.fullComm().NumProc(); ++i) {
        //  if (i == m_env.fullRank()) {
        //    std::cout << " In uqScalarFunctionSynchronizerClass<V,M>::callFunction(), just before double Bcast()"
        //              << ": fullRank "       << m_env.fullRank()
        //              << ", subEnvironment " << m_env.subId()
        //              << ", subRank "        << m_env.subRank()
        //              << ": buffer related to first double Bcast() is ready to be broadcasted"
        //              << " and has size "      << bufferDouble.size()
        //              << std::endl;
        //    if (m_env.subRank() == 0) {
	//      std::cout << "Buffer contents are";
        //      for (unsigned int i = 0; i < bufferDouble.size(); ++i) {
	//        std::cout << " " << bufferDouble[i];
        //    }
	//      std::cout << std::endl;
        //    }
        //  }
        //  m_env.fullComm().Barrier();
        //}
        //if (m_env.fullRank() == 0) std::cout << "Sleeping 3 seconds..."
        //                                 << std::endl;
        //sleep(3);

        count = (int) bufferDouble.size();
        mpiRC = MPI_Bcast ((void *) &bufferDouble[0], count, MPI_DOUBLE, 0, m_env.subComm().Comm());
        UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                            m_env.fullRank(),
                            "uqScalarFunctionSynchronizerClass<V,M>::callFunction()",
                            "failed broadcast 2 of 3");

        if (m_env.subRank() != 0) {
          V tmpVec(m_auxVec);
          for (unsigned int i = 0; i < tmpVec.sizeLocal(); ++i) {
            tmpVec[i] = bufferDouble[i];
          }
          internalValues = new V(tmpVec);
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
          mpiRC = MPI_Bcast ((void *) &bufferDouble[0], count, MPI_DOUBLE, 0, m_env.subComm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              m_env.fullRank(),
                              "uqScalarFunctionSynchronizerClass<V,M>::callFunction()",
                              "failed broadcast 3 of 3");

          if (m_env.subRank() != 0) {
            V tmpVec(m_auxVec);
            for (unsigned int i = 0; i < tmpVec.sizeLocal(); ++i) {
              tmpVec[i] = bufferDouble[i];
            }
            internalDirection = new V(tmpVec);
          }
        }

        ///////////////////////////////////////////////
        // All processors now call 'scalarFunction()'
        ///////////////////////////////////////////////
        if (m_env.subRank() != 0) {
          if (bufferChar[2] == '1') internalGrad    = new V(m_auxVec);
          if (bufferChar[3] == '1') internalHessian = new M(m_auxVec);
          if (bufferChar[4] == '1') internalEffect  = new V(m_auxVec);
        }

        m_env.syncPrintDebugMsg("In uqScalarFunctionSynchronizerClass<V,M>::callFunction(), just before actual minus2LnValue()",3,3000000,m_env.subComm());
        m_env.subComm().Barrier();
        result = m_scalarFunction.minus2LnValue(*internalValues,
                                                internalDirection,
                                                internalGrad,
                                                internalHessian,
                                                internalEffect);
      } // if (bufferChar[0] == '1')

      /////////////////////////////////////////////////
      // Prepare to exit routine or to stay in it
      /////////////////////////////////////////////////
      if (m_env.subRank() == 0) {
        stayInRoutine = false; // Always for processor 0
      }
      else {
        if (internalValues    != NULL) delete internalValues;
        if (internalDirection != NULL) delete internalDirection;
        if (internalGrad      != NULL) delete internalGrad;
        if (internalHessian   != NULL) delete internalHessian;
        if (internalEffect    != NULL) delete internalEffect;

        stayInRoutine = (vecValues == NULL) && (bufferChar[0] == '1');
        //if (!stayInRoutine) std::cout << "Fullrank " << m_env.fullRank() << " is leaving scalarFunctionSync()" << std::endl;
      }
    } while (stayInRoutine);
  }
  else {
    UQ_FATAL_TEST_MACRO(vecValues == NULL,
                        m_env.fullRank(),
                        "uqScalarFunctionSynchronizerClass<V,M>::callFunction()",
                        "vecValues should not be NULL");

    m_env.subComm().Barrier();
    result = m_scalarFunction.minus2LnValue(*vecValues,
                                            vecDirection,
                                            gradVector,
                                            hessianMatrix,
                                            hessianEffect);
  }

  return result;
}

#endif // __UQ_SCALAR_FUNCTION_SYNCHRONIZER_H__
