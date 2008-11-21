/* uq/examples/queso/pyramid/uqRadEx.h
 *
 * Copyright (C) 2008 The QUESO Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_RAD_VALIDATION_H__
#define __UQ_RAD_VALIDATION_H__

#include <uqModelValidation.h>
#include <uqVectorSubset.h>
#include <uqAsciiTable.h>

//********************************************************
// Likelihood function object for both inverse problems of the validation cycle.
// A likelihood function object is provided by user and is called by the UQ library.
// This likelihood function object consists of data and routine.
//********************************************************

// The (user defined) data class that carries the data needed by the (user defined) likelihood routine
template<class P_V, class P_M>
struct
radLikelihoodRoutine_DataClass
{
  radLikelihoodRoutine_DataClass(const uqBaseEnvironmentClass& env,
                                 const char* inpName);
 ~radLikelihoodRoutine_DataClass();

  std::vector<double> m_P;   // pressure
  std::vector<double> m_T;   // temperatures
  std::vector<double> m_I;   // integrated intensities of radiation
  std::vector<double> m_Var; //
};

template<class P_V, class P_M>
radLikelihoodRoutine_DataClass<P_V,P_M>::radLikelihoodRoutine_DataClass(
  const uqBaseEnvironmentClass& env,
  const char* inpName)
  :
  m_P(1,0.),
  m_T(1,0.),
  m_I(1,0.),
  m_Var(1,0.)
{
  // Read experimental data
  if (inpName) {
    if (env.rank() == 0) {
      std::cout << "In radLikelihoodRoutine_DataClass(), reading file '"
                << inpName << "'\n"
                << std::endl;
    }

    // Open input file on experimental data
    FILE *inp;
    inp = fopen(inpName,"r");

    // Read data
    unsigned int numLines;
    fscanf(inp,"%d",&numLines);
    m_P.resize(numLines,0.);
    m_T.resize(numLines,0.);
    m_I.resize(numLines,0.);
    m_Var.resize(numLines,0.);
    for (unsigned int i = 0; i < numLines; ++i) {
      fscanf(inp,"%lf %lf %lf %lf",&m_P[i],&m_T[i],&m_I[i],&m_Var[i]);
    }

    for (unsigned int i = 0; i < numLines; ++i) {
      printf("i = %d, m_P[i] = %lf, m_T[i] = %lf, m_I[i] = %lf, m_Var[i] = %lf\n",i,m_P[i],m_T[i],m_I[i],m_Var[i]);
    }

    // Close input file on experimental data
    fclose(inp);
  }
}

template<class P_V, class P_M>
radLikelihoodRoutine_DataClass<P_V,P_M>::~radLikelihoodRoutine_DataClass()
{
}

// The actual (user defined) likelihood routine
template<class P_V,class P_M>
double
radLikelihoodRoutine(const P_V& paramValues, const void* functionDataPtr)
{
  //DOUBLE PRECISION resultValue = 0.;

  //DOUBLE PRECISION cp1 = paramValues[0];
  //DOUBLE PRECISION cp2 = paramValues[1];
  //DOUBLE PRECISION ct1 = paramValues[2];
  //DOUBLE PRECISION ct2 = paramValues[3];
  //DOUBLE PRECISION ct3 = paramValues[4];
  //CALL SPECTRAL_MODEL(cp1,cp2,ct1,ct2,ct3,resultValue);

  double resultValue = 0.;

  double CP1 = paramValues[0];
  double CP2 = paramValues[1];
  double CT1 = paramValues[2];
  double CT2 = paramValues[3];
  double CT3 = paramValues[4];

  const std::vector<double>& m_P   = ((radLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_P;
  const std::vector<double>& m_T   = ((radLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_T;
  const std::vector<double>& m_I   = ((radLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_I;
  const std::vector<double>& m_Var = ((radLikelihoodRoutine_DataClass<P_V,P_M> *) functionDataPtr)->m_Var;

  double DS=.1016;
  double SIG=5.67e-8;

  for (unsigned int i = 0; i < m_P.size(); ++i) {
    double CP=CP1*m_P[i]-CP2;

    double K_G=CP*(CT1*m_T[i]*m_T[i]+CT2*m_T[i]+CT3);

    double IB=SIG*(m_T[i]*m_T[i]*m_T[i]*m_T[i])/M_PI;

    double I=IB*(1.-exp(-K_G*DS));

    resultValue += (I-m_I[i])*(I-m_I[i])/m_Var[i];
  }

  return resultValue;
}

// The (user defined) grad likelihood routine
template<class P_V,class P_M>
void
radLikelihoodGradRoutine(const P_V& paramValues, const void* functionDataPtr, P_V& gradVector)
{
  gradVector *= 0.;

  return;
}

// The (user defined) grad likelihood routine
template<class P_V,class P_M>
void
radLikelihoodHessianRoutine(const P_V& paramValues, const void* functionDataPtr, P_M& hessianMatrix)
{
  hessianMatrix *= 0.;

  return;
}

//********************************************************
// QoI function object for both forward problems of the validation cycle.
// A QoI function object is provided by user and is called by the UQ library.
// This QoI function object consists of data and routine.
//********************************************************
// The (user defined) data class that carries the data needed by the (user defined) qoi routine
template<class P_V,class P_M,class Q_V, class Q_M>
struct
radQoiRoutine_DataClass
{
  double m_P;
  double m_T;
};

// The actual (user defined) qoi routine
template<class P_V,class P_M,class Q_V,class Q_M>
void radQoiRoutine(const P_V& paramValues, const void* functionDataPtr, Q_V& qoiValues)
{
  double CP1 = paramValues[0];
  double CP2 = paramValues[1];
  double CT1 = paramValues[2];
  double CT2 = paramValues[3];
  double CT3 = paramValues[4];

  double m_P  = ((radQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_P;
  double m_T  = ((radQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M> *) functionDataPtr)->m_T;

  double DS=.1016;
  double SIG=5.67e-8;

    double CP=CP1*m_P-CP2;

    double K_G=CP*(CT1*m_T*m_T+CT2*m_T+CT3);

    double IB=SIG*(m_T*m_T*m_T*m_T)/M_PI;

    double I=IB*(1.-exp(-K_G*DS));

    qoiValues[0] = M_PI*I;
	
  //if ((paramValues.env().verbosity() >= 3) && (paramValues.env().rank() == 0)) {
  //  printf("In radQoiRoutine(), A = %g, E = %g, beta = %.3lf, criticalTime = %.3lf, criticalMass = %.3lf: qoi = %lf.\n",A,E,beta,criticalTime,criticalMass,qoiValues[0]);
  //}

  return;
}

//********************************************************
// Class "uqRadValidation", instantiated by main()
//********************************************************
template <class P_V,class P_M,class Q_V,class Q_M>
class uqRadValidationClass : public uqModelValidationClass<P_V,P_M,Q_V,Q_M>
{
public:
  uqRadValidationClass(const uqBaseEnvironmentClass& env,
                       const char*               prefix);
 ~uqRadValidationClass();

  void run();

private:
  void  runCalibrationStage();
  void  runValidationStage();
  void  runComparisonStage();

  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_env;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_prefix;
  using uqModelValidationClass<P_V,P_M,Q_V,Q_M>::m_cycle;

  uqAsciiTableClass<P_V,P_M>*               m_paramsTable;
  const EpetraExt::DistArray<std::string>*  m_paramNames;         // instantiated outside this class!!
  P_V*                                      m_paramMinValues;     // instantiated outside this class!!
  P_V*                                      m_paramMaxValues;     // instantiated outside this class!!
  P_V*                                      m_paramInitialValues; // instantiated outside this class!!
  uqVectorSpaceClass<P_V,P_M>*              m_paramSpace;
  uqVectorSetClass<P_V,P_M>*                m_paramDomain;

  uqAsciiTableClass<P_V,P_M>*               m_qoisTable;
  const EpetraExt::DistArray<std::string>*  m_qoiNames; // instantiated outside this class!!
  uqVectorSpaceClass<Q_V,Q_M>*              m_qoiSpace;

  double                                    m_predT;
  double                                    m_predP;

  uqBaseVectorRVClass<P_V,P_M>*             m_calPriorRv;
  radLikelihoodRoutine_DataClass<P_V,P_M>*  m_calLikelihoodRoutine_Data;
  uqBaseScalarFunctionClass<P_V,P_M>*       m_calLikelihoodFunctionObj;
  radQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>* m_calQoiRoutine_Data;

  radLikelihoodRoutine_DataClass<P_V,P_M>*  m_valLikelihoodRoutine_Data;
  uqBaseScalarFunctionClass<P_V,P_M>*       m_valLikelihoodFunctionObj;
  radQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>* m_valQoiRoutine_Data;
};

template <class P_V,class P_M,class Q_V,class Q_M>
uqRadValidationClass<P_V,P_M,Q_V,Q_M>::uqRadValidationClass(
  const uqBaseEnvironmentClass& env,
  const char*                   prefix)
  :
  uqModelValidationClass<P_V,P_M,Q_V,Q_M>(env,prefix),
  m_paramsTable              (NULL),
  m_paramNames               (NULL),
  m_paramMinValues           (NULL),
  m_paramMaxValues           (NULL),
  m_paramInitialValues       (NULL),
  m_paramSpace               (NULL),
  m_paramDomain              (NULL),
  m_qoisTable                (NULL),
  m_qoiNames                 (NULL),
  m_qoiSpace                 (NULL),
  m_calPriorRv               (NULL),
  m_calLikelihoodRoutine_Data(NULL),
  m_calLikelihoodFunctionObj (NULL),
  m_calQoiRoutine_Data       (NULL),
  m_valLikelihoodRoutine_Data(NULL),
  m_valLikelihoodFunctionObj (NULL),
  m_valQoiRoutine_Data       (NULL)
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqRadValidation::constructor()\n"
              << std::endl;
  }

  // Read Ascii file with information on parameters
  m_paramsTable = new uqAsciiTableClass<P_V,P_M> (m_env,
                                                  5,    // # of rows
                                                  3,    // # of cols after 'parameter name': min + max + initial value for Markov chain
                                                  NULL, // All extra columns are of 'double' type
                                                  "params.tab");

  m_paramNames = &(m_paramsTable->stringColumn(0));
  m_paramMinValues     = new P_V(m_paramsTable->doubleColumn(1));
  m_paramMaxValues     = new P_V(m_paramsTable->doubleColumn(2));
  m_paramInitialValues = new P_V(m_paramsTable->doubleColumn(3));

  m_paramSpace = new uqVectorSpaceClass<P_V,P_M>(m_env,
                                                 "param_", // Extra prefix before the default "space_" prefix
                                                 m_paramsTable->numRows(),
                                                 m_paramNames);

  m_paramDomain = new uqBoxSubsetClass<P_V,P_M>("param_",
                                                *m_paramSpace,
                                                *m_paramMinValues,
                                                *m_paramMaxValues);

  // Read Ascii file with information on qois
  m_qoisTable = new uqAsciiTableClass<P_V,P_M>(m_env,
                                               1,    // # of rows
                                               0,    // # of cols after 'parameter name': none
                                               NULL, // All extra columns are of 'double' type
                                               "qois.tab");

  m_qoiNames = &(m_qoisTable->stringColumn(0));

  m_qoiSpace = new uqVectorSpaceClass<Q_V,Q_M>(m_env,
                                               "qoi_", // Extra prefix before the default "space_" prefix
                                               m_qoisTable->numRows(),
                                               m_qoiNames);

  // Instantiate the validation cycle
  m_cycle = new uqValidationCycleClass<P_V,P_M,Q_V,Q_M>(m_env,
                                                        m_prefix.c_str(), // Use the prefix passed above
                                                        *m_paramSpace,
                                                        *m_qoiSpace);

  m_predT = 10337.;
  m_predP = 0.77;

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqRadValidation::constructor()\n"
              << std::endl;
  }

  return;
}

template <class P_V,class P_M,class Q_V,class Q_M>
uqRadValidationClass<P_V,P_M,Q_V,Q_M>::~uqRadValidationClass()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqRadValidation::destructor()"
              << std::endl;
  }

  if (m_valQoiRoutine_Data)        delete m_valQoiRoutine_Data;
  if (m_valLikelihoodFunctionObj)  delete m_valLikelihoodFunctionObj;
  if (m_valLikelihoodRoutine_Data) delete m_valLikelihoodRoutine_Data;
  if (m_calQoiRoutine_Data)        delete m_calQoiRoutine_Data;
  if (m_calLikelihoodFunctionObj)  delete m_calLikelihoodFunctionObj;
  if (m_calLikelihoodRoutine_Data) delete m_calLikelihoodRoutine_Data;
  if (m_calPriorRv)                delete m_calPriorRv;
  if (m_qoiSpace)                  delete m_qoiSpace;
//if (m_qoiNames)                  delete m_qoiNames; // instantiated outside this class!!
  if (m_qoisTable)                 delete m_qoisTable;
  if (m_paramDomain)               delete m_paramDomain;
  if (m_paramSpace)                delete m_paramSpace;
//if (m_paramInitialValues)        delete m_paramInitialValues; // instantiated outside this class!!
//if (m_paramMaxValues)            delete m_paramMaxValues;     // instantiated outside this class!!
//if (m_paramMinValues)            delete m_paramMinValues;     // instantiated outside this class!!
//if (m_paramNames)                delete m_paramNames;         // instantiated outside this class!!
  if (m_paramsTable)               delete m_paramsTable;

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqRadValidation::destructor()"
              << std::endl;
  }
}

template <class P_V,class P_M,class Q_V,class Q_M>
void
uqRadValidationClass<P_V,P_M,Q_V,Q_M>::run()
{
  if (m_env.rank() == 0) {
    std::cout << "Entering uqRadValidation::run()"
              << std::endl;
  }

  runCalibrationStage();
  runValidationStage();
  runComparisonStage();

  if (m_env.rank() == 0) {
    std::cout << "Leaving uqRadValidation::run()"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqRadValidationClass<P_V,P_M,Q_V,Q_M>::runCalibrationStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqRadValidation::runCalibrationStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  m_calPriorRv = new uqUniformVectorRVClass<P_V,P_M> ("cal_prior_", // Extra prefix before the default "rv_" prefix
                                                      *m_paramDomain);

  m_calLikelihoodRoutine_Data = new radLikelihoodRoutine_DataClass<P_V,P_M>(m_env,
                                                                            "I1.dat");
  m_calLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("cal_like_",
                                                                         *m_paramDomain,
                                                                         radLikelihoodRoutine<P_V,P_M>,
                                                                         radLikelihoodGradRoutine<P_V,P_M>,
                                                                         radLikelihoodHessianRoutine<P_V,P_M>,
                                                                         (void *) m_calLikelihoodRoutine_Data,
                                                                         true); // the routine computes [-2.*ln(function)]

  m_cycle->setCalIP(*m_calPriorRv,
                    *m_calLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* calProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(m_cycle->calIP().priorRv().pdf().domainVarVector(),
                                                                                                  *m_paramInitialValues);
  m_cycle->calIP().solveWithBayesMarkovChain(*m_paramInitialValues,
                                             calProposalCovMatrix);
  delete calProposalCovMatrix;

  // Deal with forward problem
  m_calQoiRoutine_Data = new radQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>();
  m_calQoiRoutine_Data->m_T = m_predT;
  m_calQoiRoutine_Data->m_P = m_predP;

  m_cycle->setCalFP(radQoiRoutine<P_V,P_M,Q_V,Q_M>,
                    (void *) m_calQoiRoutine_Data);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  m_cycle->calFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqRadValidation::runCalibrationStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqRadValidation::runCalibrationStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqRadValidationClass<P_V,P_M,Q_V,Q_M>::runValidationStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqRadValidation::runValidationStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  // Deal with inverse problem
  m_valLikelihoodRoutine_Data = new radLikelihoodRoutine_DataClass<P_V,P_M>(m_env,
                                                                            "I2.dat");

  m_valLikelihoodFunctionObj = new uqGenericScalarFunctionClass<P_V,P_M>("val_like_",
                                                                         *m_paramDomain,
                                                                         radLikelihoodRoutine<P_V,P_M>,
                                                                         radLikelihoodGradRoutine<P_V,P_M>,
                                                                         radLikelihoodHessianRoutine<P_V,P_M>,
                                                                         (void *) m_valLikelihoodRoutine_Data,
                                                                         true); // the routine computes [-2.*ln(function)]

  m_cycle->setValIP(*m_valLikelihoodFunctionObj);

  // Solve inverse problem = set 'pdf' and 'realizer' of 'postRv'
  P_M* valProposalCovMatrix = m_cycle->calIP().postRv().imageSet().vectorSpace().newGaussianMatrix(m_cycle->calIP().postRv().realizer().imageVarVector(),  // Use 'realizer()' because the posterior rv was computed with Markov Chain
                                                                                                   m_cycle->calIP().postRv().realizer().imageExpVector()); // Use these values as the initial values
  m_cycle->valIP().solveWithBayesMarkovChain(m_cycle->calIP().postRv().realizer().imageExpVector(),
                                             valProposalCovMatrix);
  delete valProposalCovMatrix;

  // Deal with forward problem
  m_valQoiRoutine_Data = new radQoiRoutine_DataClass<P_V,P_M,Q_V,Q_M>();
  m_valQoiRoutine_Data->m_T = m_predT;
  m_valQoiRoutine_Data->m_P = m_predP;

  m_cycle->setValFP(radQoiRoutine<P_V,P_M,Q_V,Q_M>,
                    (void *) m_valQoiRoutine_Data);

  // Solve forward problem = set 'realizer' and 'cdf' of 'qoiRv'
  m_cycle->valFP().solveWithMonteCarlo(); // no extra user entities needed for Monte Carlo algorithm

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqRadValidation::runValidationStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqRadValidation::runValidationStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}

template<class P_V,class P_M,class Q_V,class Q_M>
void 
uqRadValidationClass<P_V,P_M,Q_V,Q_M>::runComparisonStage()
{
  int iRC;
  struct timeval timevalRef;
  struct timeval timevalNow;

  iRC = gettimeofday(&timevalRef, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Entering uqRadValidation::runComparisonStage() at " << ctime(&timevalRef.tv_sec)
              << std::endl;
  }

  if (m_cycle->calFP().computeSolutionFlag() &&
      m_cycle->valFP().computeSolutionFlag()) {
    Q_V* epsilonVec = m_cycle->calFP().qoiRv().imageSet().vectorSpace().newVector(0.02);
    Q_V cdfDistancesVec(m_cycle->calFP().qoiRv().imageSet().vectorSpace().zeroVector());
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Test independence of 'distance' w.r.t. order of cdfs
    horizontalDistances(m_cycle->valFP().qoiRv().cdf(),
                        m_cycle->calFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "                             << *epsilonVec
                << ", cdfDistancesVec (switched order of cdfs) = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.04
    epsilonVec->cwSet(0.04);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.06
    epsilonVec->cwSet(0.06);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.08
    epsilonVec->cwSet(0.08);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    // Epsilon = 0.10
    epsilonVec->cwSet(0.10);
    horizontalDistances(m_cycle->calFP().qoiRv().cdf(),
                        m_cycle->valFP().qoiRv().cdf(),
                        *epsilonVec,
                        cdfDistancesVec);
    if (m_cycle->env().rank() == 0) {
      std::cout << "For epsilonVec = "    << *epsilonVec
                << ", cdfDistancesVec = " << cdfDistancesVec
                << std::endl;
    }

    delete epsilonVec;
  }

  iRC = gettimeofday(&timevalNow, NULL);
  if (m_env.rank() == 0) {
    std::cout << "Leaving uqRadValidation::runComparisonStage() at " << ctime(&timevalNow.tv_sec)
              << "Total uqRadValidation::runComparisonStage() run time = " << timevalNow.tv_sec - timevalRef.tv_sec
              << " seconds"
              << std::endl;
  }

  return;
}
#endif // __UQ_RAD_VALIDATION_H__
