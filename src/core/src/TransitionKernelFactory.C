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

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/TransitionKernelFactory.h>
#include <queso/ScaledCovMatrixTKGroup.h>
#include <queso/TransformedScaledCovMatrixTKGroup.h>
#include <queso/HessianCovMatricesTKGroup.h>

namespace QUESO
{

template <>
std::map<std::string, Factory<BaseTKGroup<GslVector, GslMatrix> > *> &
Factory<BaseTKGroup<GslVector, GslMatrix> >::factory_map()
{
  static std::map<std::string, Factory<BaseTKGroup<GslVector, GslMatrix> > *> _factory_map;

  return _factory_map;
}

// SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type
// TransitionKernelFactory::build_tk(
//     const MhOptionsValues & options,
//     const VectorSpace<GslVector, GslMatrix> & v,
//     const std::vector<double> & dr_scales,
//     const ScalarFunctionSynchronizer<GslVector, GslMatrix> & pdf_synchronizer,
//     GslMatrix & initial_cov_matrix,
//     const BaseJointPdf<GslVector, GslMatrix> & target_pdf)
// {
//   SharedPtr<BaseTKGroup<GslVector, GslMatrix> >::Type new_tk;
//
//   if (options.m_tkUseLocalHessian) {
//     new_tk.reset(new HessianCovMatricesTKGroup<GslVector,GslMatrix>(
//           options.m_prefix.c_str(),
//           v,
//           dr_scales,
//           pdf_synchronizer));
//   }
//   else {
//     if (options.m_initialProposalCovMatrixDataInputFileName != ".") {
//       std::set<unsigned int> tmpSet;
//       tmpSet.insert(v.env().subId());
//       initial_cov_matrix.subReadContents((options.m_initialProposalCovMatrixDataInputFileName+"_sub"+v.env().subIdString()),
//                                           options.m_initialProposalCovMatrixDataInputFileType,
//                                           tmpSet);
//     }
//
//     // Decide whether or not to do logit transform
//     if (options.m_doLogitTransform) {
//       // Variable transform of initial proposal cov matrix has already been
//       // done
//
//       // We need this dynamic_cast to BoxSubset so that m_tk can inspect the
//       // domain bounds and do the necessary transform
//       new_tk.reset(new TransformedScaledCovMatrixTKGroup<GslVector, GslMatrix>(
//           options.m_prefix.c_str(),
//           dynamic_cast<const BoxSubset<GslVector, GslMatrix> & >(target_pdf.domainSet()),
//           dr_scales, initial_cov_matrix));
//     }
//     else {
//       new_tk.reset(new ScaledCovMatrixTKGroup<GslVector, GslMatrix>(
//           options.m_prefix.c_str(), v, dr_scales,
//           initial_cov_matrix));
//     }
//   }
//
//   return new_tk;
// }

const VectorSpace<GslVector, GslMatrix> * TransitionKernelFactory::m_vectorSpace = NULL;

const std::vector<double> * TransitionKernelFactory::m_dr_scales = NULL;

const ScalarFunctionSynchronizer<GslVector, GslMatrix> * TransitionKernelFactory::m_pdf_synchronizer = NULL;

GslMatrix * TransitionKernelFactory::m_initial_cov_matrix = NULL;

const MhOptionsValues * TransitionKernelFactory::m_options = NULL;

const BaseJointPdf<GslVector, GslMatrix> * TransitionKernelFactory::m_target_pdf = NULL;

} // namespace QUESO
