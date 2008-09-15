/* uq/libs/mcmc/inc/uqProblem.h
 *
 * Copyright (C) 2008 The PECOS Team, http://queso.ices.utexas.edu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __UQ_PROBLEM_H__
#define __UQ_PROBLEM_H__

#include <uqProblemSlice.h>

#define UQ_PROBLEM_PREFIX_FOR_NO_SLICE_PREFIX "."

// _ODV = option default value
#define UQ_PROBLEM_NUM_SLICES_ODV     0
#define UQ_PROBLEM_SLICE_PREFIXES_ODV UQ_PROBLEM_PREFIX_FOR_NO_SLICE_PREFIX
#define UQ_PROBLEM_SLICE_ORDER_ODV    "0"

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
class uqProblemClass
{
public:
  uqProblemClass(const uqEnvironmentClass& env);
 ~uqProblemClass();

        void                                          instantiateSlice (unsigned int                                           sliceId,
                                                                        const uqProbDensity_BaseClass       <P_V,P_M>*         m2lPriorParamDensityObj,  // Set in substep x.1 in applications setting a problem slice
                                                                        const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>* m2lLikelihoodFunctionObj, // Set in substep x.2
                                                                        P_M*                                                   proposalCovMatrix,        // Set in substep x.3
                                                                        const uqProposalDensity_BaseClass   <P_V,P_M>*         proposalDensityObj,       // Set in substep x.3
                                                                        const uqProposalGenerator_BaseClass <P_V,P_M>*         proposalGeneratorObj,     // Set in substep x.3
                                                                        const uqProbDensity_BaseClass       <P_V,P_M>*         propagParamDensityObj,    // Set in substep x.4
                                                                        const uqSampleGenerator_BaseClass   <P_V,P_M>*         propagParamGeneratorObj,  // Set in substep x.4
                                                                        const uqQoIFunction_BaseClass       <P_V,P_M,Q_V,Q_M>* qoiFunctionObj);          // Set in substep x.5

  const uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& slice            (unsigned int sliceId) const;
        void                                          quantify         ();

        void                                          print            (std::ostream& os) const;

private:
        void                                          defineMyOptions  (po::options_description& optionsDesc);
        void                                          getMyOptionValues(po::options_description& optionsDesc);

  const uqEnvironmentClass&                                  m_env;

  po::options_description*                                   m_optionsDesc;
  std::string                                                m_option_help;
  std::string                                                m_option_numSlices;
  std::string                                                m_option_slicePrefixes;
  std::string                                                m_option_sliceOrder;

  unsigned int                                               m_numSlices;
  std::vector<std::string>                                   m_slicePrefixes;
  std::vector<unsigned int>                                  m_sliceOrder;
  std::vector<uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>*> m_slices;
};

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj);

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::uqProblemClass(const uqEnvironmentClass& env)
  :
  m_env                 (env),
  m_optionsDesc         (new po::options_description("Uncertainty Quantification Problem")),
  m_option_help         ("uqProblem_help"         ),
  m_option_numSlices    ("uqProblem_numSlices"    ),
  m_option_slicePrefixes("uqProblem_slicePrefixes"),
  m_option_sliceOrder   ("uqProblem_sliceOrder"   ),
  m_numSlices           (UQ_PROBLEM_NUM_SLICES_ODV),
  m_slicePrefixes       (1,UQ_PROBLEM_SLICE_PREFIXES_ODV),
  m_sliceOrder          (0),
  m_slices              (0)
{
  if (m_env.rank() == 0) std::cout << "Entering uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << std::endl;

  defineMyOptions                (*m_optionsDesc);
  m_env.scanInputFileForMyOptions(*m_optionsDesc);
  getMyOptionValues              (*m_optionsDesc);

  if (m_env.rank() == 0) std::cout << "In uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << ": after getting values of options, state of object is:"
                                   << "\n" << *this
                                   << std::endl;

  m_slices.resize(m_numSlices,NULL);

  if (m_env.rank() == 0) std::cout << "Leaving uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::constructor()"
                                   << std::endl;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::~uqProblemClass()
{
  for (unsigned int i = 0; i < m_slices.size(); ++i) {
    if (m_slices[i] != NULL) delete m_slices[i];
  }
  if (m_optionsDesc) delete m_optionsDesc;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::defineMyOptions(
  po::options_description& optionsDesc)
{
  optionsDesc.add_options()
    (m_option_help.c_str(),                                                                                   "produce help message for UQ problem"    )
    (m_option_numSlices.c_str(),     po::value<unsigned int>()->default_value(UQ_PROBLEM_NUM_SLICES_ODV),     "number of slices"                       )
    (m_option_slicePrefixes.c_str(), po::value<std::string >()->default_value(UQ_PROBLEM_SLICE_PREFIXES_ODV), "list of prefix(es) for slice(s)"        )
    (m_option_sliceOrder.c_str(),    po::value<std::string >()->default_value(UQ_PROBLEM_SLICE_ORDER_ODV),    "order of slices during solution process")
  ;

  return;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
  uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::getMyOptionValues(
  po::options_description& optionsDesc)
{
  if (m_env.allOptionsMap().count(m_option_help.c_str())) {
    std::cout << optionsDesc
              << std::endl;
  }

  if (m_env.allOptionsMap().count(m_option_numSlices.c_str())) {
    m_numSlices = m_env.allOptionsMap()[m_option_numSlices.c_str()].as<unsigned int>();
  }

  if (m_env.allOptionsMap().count(m_option_slicePrefixes.c_str())) {
    m_slicePrefixes.clear();
    std::string inputString = m_env.allOptionsMap()[m_option_slicePrefixes.c_str()].as<std::string>();
    uqMiscReadWordsFromString(inputString,m_slicePrefixes);
    //std::cout << "In uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::getMyOptionValues()"
    //          << ": m_slicePrefixes[0] = " << m_slicePrefixes[0]
    //          << std::endl;
  }

  UQ_FATAL_TEST_MACRO(m_numSlices != m_slicePrefixes.size(),
                      m_env.rank(),
                      "uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::getMyOptionsValues()",
                      "option 'numSlices' is not equal to size of array of slice prefixes");

  if (m_env.allOptionsMap().count(m_option_sliceOrder.c_str())) {
    m_sliceOrder.clear();
    std::vector<double> tmpOrder(0,0.);
    std::string inputString = m_env.allOptionsMap()[m_option_sliceOrder.c_str()].as<std::string>();
    uqMiscReadDoublesFromString(inputString,tmpOrder);

    if (tmpOrder.size() > 0) {
      m_sliceOrder.resize(tmpOrder.size(),0);
      for (unsigned int i = 0; i < m_sliceOrder.size(); ++i) {
        m_sliceOrder[i] = (unsigned int) tmpOrder[i];
      }
    }
  }
  else {
    m_sliceOrder.resize(m_numSlices,0);
    for (unsigned int i = 0; i < m_sliceOrder.size(); ++i) {
      m_sliceOrder[i] = i;
    }
  }

  UQ_FATAL_TEST_MACRO(m_numSlices != m_sliceOrder.size(),
                      m_env.rank(),
                      "uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::getMyOptionsValues()",
                      "option 'numSlices' is not equal to size of array of slice order");

  for (unsigned int i = 0; i < m_sliceOrder.size(); ++i) {
    UQ_FATAL_TEST_MACRO(m_sliceOrder[i] >= m_numSlices,
                        m_env.rank(),
                        "uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::getMyOptionsValues()",
                        "a slice id on the array of slice order is too big");
  }

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::instantiateSlice(
  unsigned int                                           sliceId,
  const uqProbDensity_BaseClass       <P_V,P_M>*         m2lPriorParamDensityObj,  // Set in substep x.1
  const uqLikelihoodFunction_BaseClass<P_V,P_M,L_V,L_M>* m2lLikelihoodFunctionObj, // Set in substep x.2
  P_M*                                                   proposalCovMatrix,        // Set in substep x.3
  const uqProposalDensity_BaseClass   <P_V,P_M>*         proposalDensityObj,       // Set in substep x.3
  const uqProposalGenerator_BaseClass <P_V,P_M>*         proposalGeneratorObj,     // Set in substep x.3
  const uqProbDensity_BaseClass       <P_V,P_M>*         propagParamDensityObj,    // Set in substep x.4
  const uqSampleGenerator_BaseClass   <P_V,P_M>*         propagParamGeneratorObj,  // Set in substep x.4
  const uqQoIFunction_BaseClass       <P_V,P_M,Q_V,Q_M>* qoiFunctionObj)           // Set in substep x.5
{
  m_slices[sliceId] = new uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>(m_env,
                                                                       m_slicePrefixes[sliceId].c_str(),
                                                                       m2lPriorParamDensityObj,
                                                                       m2lLikelihoodFunctionObj,
                                                                       proposalCovMatrix,
                                                                       proposalDensityObj,
                                                                       proposalGeneratorObj,
                                                                       propagParamDensityObj,
                                                                       propagParamGeneratorObj,
                                                                       qoiFunctionObj);

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
const uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>&
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::slice(unsigned int sliceId) const
{
  return *(m_slices[sliceId]);
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::quantify()
{
  for (unsigned int i = 0; i < m_numSlices; ++i) {
    uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>* currentSlice = m_slices[m_sliceOrder[i]];
    if (currentSlice->isCalibRequested()) {
      if (currentSlice->calibInputSliceId() < 0) {
        currentSlice->calibrate();
      }
      else {
        unsigned int inputSliceId = (unsigned int) currentSlice->calibInputSliceId();
        uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>* inputSlice = m_slices[inputSliceId];
        currentSlice->calibrate(inputSlice->posteriorParamDensityObj());
      }
    }
  }

  for (unsigned int i = 0; i < m_numSlices; ++i) {
    uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>* currentSlice = m_slices[m_sliceOrder[i]];
    if (currentSlice->isPropagRequested()) {
      if (currentSlice->propagInputSliceId() < 0) {
        currentSlice->propagate();
      }
      else {
        unsigned int inputSliceId = (unsigned int) currentSlice->propagInputSliceId();
        uqProblemSliceClass<P_V,P_M,L_V,L_M,Q_V,Q_M>* inputSlice = m_slices[inputSliceId];
        currentSlice->propagate(inputSlice->posteriorParamGeneratorObj());
      }
    }
  }

  return;
}

template <class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
void
uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>::print(std::ostream& os) const
{
  os << "\n" << m_option_numSlices     << " = " << m_numSlices
     << "\n" << m_option_slicePrefixes << " = ";
  for (unsigned int i = 0; i < m_slicePrefixes.size(); ++i) {
    os << m_slicePrefixes[i] << " ";
  }
  os << "\n" << m_option_sliceOrder << " = ";
  for (unsigned int i = 0; i < m_sliceOrder.size(); ++i) {
    os << m_sliceOrder[i] << " ";
  }
  os << std::endl;
}

template<class P_V,class P_M,class L_V,class L_M,class Q_V,class Q_M>
std::ostream& operator<<(std::ostream& os, const uqProblemClass<P_V,P_M,L_V,L_M,Q_V,Q_M>& obj)
{
  obj.print(os);

  return os;
}
#endif // __UQ_PROBLEM_H__
