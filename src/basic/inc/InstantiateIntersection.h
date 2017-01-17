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

#ifndef UQ_INSTANTIATE_INTERSECTION_H
#define UQ_INSTANTIATE_INTERSECTION_H

#include <queso/VectorSet.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!\file InstantiateIntersection.h
 * \brief A templated method to calculate intersection of two domains (vector spaces).
 */

//! This method calculates the intersection of \c domain1 and \c domain2.
/*! It is used, for instance, to calculate the domain of the Posterior PDF, which is
 * the intersection of the domain of the Prior PDF and of the likelihood function.*/
template <class V, class M>
VectorSet<V,M>* InstantiateIntersection(const VectorSet<V,M>& domain1,
    const VectorSet<V,M>& domain2);

}  // End namespace QUESO

#endif // UQ_INSTANTIATE_INTERSECTION_H

