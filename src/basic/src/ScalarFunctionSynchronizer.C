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

#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/ScalarFunctionSynchronizer.h>
#include <queso/BayesianJointPdf.h>
#include <queso/JointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Default constructor
template <class V, class M>
ScalarFunctionSynchronizer<V,M>::ScalarFunctionSynchronizer(
    const BaseScalarFunction<V,M>& inputFunction,
    const V& auxVec)
  : m_env(inputFunction.domainSet().env()),
    m_scalarFunction(inputFunction),
    m_bayesianJointPdfPtr(dynamic_cast<const BayesianJointPdf<V,M>* >(&m_scalarFunction)),
    m_auxVec(auxVec)
{
}

// Destructor
template <class V, class M>
ScalarFunctionSynchronizer<V,M>::~ScalarFunctionSynchronizer()
{
}

// Math methods
template<class V,class M>
const VectorSet<V,M>& ScalarFunctionSynchronizer<V,M>::domainSet() const
{
  return m_scalarFunction.domainSet();
}

// Sync methods
template <class V,class M>
double ScalarFunctionSynchronizer<V,M>::callFunction(const V* vecValues,
    const V* vecDirection,
    V* gradVector,
    M* hessianMatrix,
    V* hessianEffect,
    double* extraOutput1,
    double* extraOutput2) const
{
  double result = 0.;

  if ((m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) &&
      (m_auxVec.numOfProcsForStorage() == 1                                  )) {
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

      m_env.subComm().syncPrintDebugMsg("In ScalarFunctionSynchronizer<V,M>::callFunction(), just before char Bcast()",3,3000000);
      //if (m_env.subId() != 0) while (true) sleep(1);

      int count = (int) bufferChar.size();
      m_env.subComm().Bcast((void *) &bufferChar[0], count, RawValue_MPI_CHAR, 0,
                            "ScalarFunctionSynchronizer<V,M>::callFunction()",
                            "failed broadcast 1 of 3");

      m_env.subComm().syncPrintDebugMsg("In ScalarFunctionSynchronizer<V,M>::callFunction(), just after char Bcast()",3,3000000);
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
        //    std::cout << " In ScalarFunctionSynchronizer<V,M>::callFunction(), just before double Bcast()"
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
        m_env.subComm().Bcast((void *) &bufferDouble[0], count, RawValue_MPI_DOUBLE, 0,
                              "ScalarFunctionSynchronizer<V,M>::callFunction()",
                              "failed broadcast 2 of 3");

        if (m_env.subRank() != 0) {
          V tmpVec(m_auxVec);
          for (unsigned int i = 0; i < tmpVec.sizeLocal(); ++i) {
            tmpVec[i] = bufferDouble[i];
          }
          internalValues = new V(tmpVec);
          //if (vecValues) *vecValues = tmpVec; // prudencio 2010-08-01
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
                                "ScalarFunctionSynchronizer<V,M>::callFunction()",
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

        m_env.subComm().syncPrintDebugMsg("In ScalarFunctionSynchronizer<V,M>::callFunction(), just before actual lnValue()",3,3000000);
        m_env.subComm().Barrier();
        result = m_scalarFunction.lnValue(*internalValues,   // input
                                          internalDirection, // input
                                          internalGrad,    // output
                                          internalHessian, // output
                                          internalEffect); // output
        if (extraOutput1) {
          if (m_bayesianJointPdfPtr) {
            *extraOutput1 = m_bayesianJointPdfPtr->lastComputedLogPrior();
          }
        }
        if (extraOutput2) {
          if (m_bayesianJointPdfPtr) {
            *extraOutput2 = m_bayesianJointPdfPtr->lastComputedLogLikelihood();
          }
        }
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
    queso_require_msg(vecValues, "vecValues should not be NULL");

    m_env.subComm().Barrier();
    result = m_scalarFunction.lnValue(*vecValues,
                                      vecDirection,
                                      gradVector,
                                      hessianMatrix,
                                      hessianEffect);
    if (extraOutput1) {
      if (m_bayesianJointPdfPtr) {
        *extraOutput1 = m_bayesianJointPdfPtr->lastComputedLogPrior();
      }
    }
    if (extraOutput2) {
      if (m_bayesianJointPdfPtr) {
        *extraOutput2 = m_bayesianJointPdfPtr->lastComputedLogLikelihood();
      }
    }
  }

  return result;
}

// Sync methods
template <class V,class M>
double ScalarFunctionSynchronizer<V,M>::callFunction(const V* vecValues,
    double* extraOutput1,
    double* extraOutput2) const
{
  double result = 0.;

  if ((m_env.numSubEnvironments() < (unsigned int) m_env.fullComm().NumProc()) &&
      (m_auxVec.numOfProcsForStorage() == 1                                  )) {
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

        if (internalValues    != NULL) bufferChar[0] = '1';
        if (internalDirection != NULL) bufferChar[1] = '1';
        if (internalGrad      != NULL) bufferChar[2] = '1';
        if (internalHessian   != NULL) bufferChar[3] = '1';
        if (internalEffect    != NULL) bufferChar[4] = '1';
      }

      m_env.subComm().syncPrintDebugMsg("In ScalarFunctionSynchronizer<V,M>::callFunction(), just before char Bcast()",3,3000000);
      //if (m_env.subId() != 0) while (true) sleep(1);

      int count = (int) bufferChar.size();
      m_env.subComm().Bcast((void *) &bufferChar[0], count, RawValue_MPI_CHAR, 0,
                            "ScalarFunctionSynchronizer<V,M>::callFunction()",
                            "failed broadcast 1 of 3");

      m_env.subComm().syncPrintDebugMsg("In ScalarFunctionSynchronizer<V,M>::callFunction(), just after char Bcast()",3,3000000);
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
        //    std::cout << " In ScalarFunctionSynchronizer<V,M>::callFunction(), just before double Bcast()"
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
        m_env.subComm().Bcast((void *) &bufferDouble[0], count, RawValue_MPI_DOUBLE, 0,
                              "ScalarFunctionSynchronizer<V,M>::callFunction()",
                              "failed broadcast 2 of 3");

        if (m_env.subRank() != 0) {
          V tmpVec(m_auxVec);
          for (unsigned int i = 0; i < tmpVec.sizeLocal(); ++i) {
            tmpVec[i] = bufferDouble[i];
          }
          internalValues = new V(tmpVec);
          //if (vecValues) *vecValues = tmpVec; // prudencio 2010-08-01
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
                                "ScalarFunctionSynchronizer<V,M>::callFunction()",
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

        m_env.subComm().syncPrintDebugMsg("In ScalarFunctionSynchronizer<V,M>::callFunction(), just before actual lnValue()",3,3000000);
        m_env.subComm().Barrier();
        result = m_scalarFunction.lnValue(*internalValues);
        if (extraOutput1) {
          if (m_bayesianJointPdfPtr) {
            *extraOutput1 = m_bayesianJointPdfPtr->lastComputedLogPrior();
          }
        }
        if (extraOutput2) {
          if (m_bayesianJointPdfPtr) {
            *extraOutput2 = m_bayesianJointPdfPtr->lastComputedLogLikelihood();
          }
        }
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
    queso_require_msg(vecValues, "vecValues should not be NULL");

    m_env.subComm().Barrier();
    result = m_scalarFunction.lnValue(*vecValues);
    if (extraOutput1) {
      if (m_bayesianJointPdfPtr) {
        *extraOutput1 = m_bayesianJointPdfPtr->lastComputedLogPrior();
      }
    }
    if (extraOutput2) {
      if (m_bayesianJointPdfPtr) {
        *extraOutput2 = m_bayesianJointPdfPtr->lastComputedLogLikelihood();
      }
    }
  }

  return result;
}

}  // End namespace QUESO

template class QUESO::ScalarFunctionSynchronizer<QUESO::GslVector, QUESO::GslMatrix>;
