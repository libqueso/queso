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

#include <queso/VectorSpace.h>
#include <queso/VectorSubset.h>
#include <queso/BoxSubset.h>
#include <queso/IntersectionSubset.h>
#include <queso/InstantiateIntersection.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

template<class V, class M>
VectorSet<V,M>*
InstantiateIntersection(const VectorSet<V,M>& domain1, const VectorSet<V,M>& domain2)
{
  VectorSet<V,M>* result = NULL;

  unsigned int dim1 = domain1.vectorSpace().dimGlobal();
  unsigned int dim2 = domain2.vectorSpace().dimGlobal();

  if (result == NULL) {
    const VectorSpace<V,M>* tmp1 = dynamic_cast<const VectorSpace<V,M>* >(&domain1);
    const VectorSpace<V,M>* tmp2 = dynamic_cast<const VectorSpace<V,M>* >(&domain2);

    if ((tmp1 != NULL) && (tmp2 != NULL)) {
      if (dim1 < dim2) {
        result = new VectorSpace<V,M>(tmp1->env(),
                                             tmp1->prefix().c_str(),
                                             tmp1->dimGlobal(),
                                             NULL);//tmp1->componentsNames());
      }
      else if (dim1 == dim2) {
        result = new VectorSpace<V,M>(tmp1->env(),
                                             tmp1->prefix().c_str(),
                                             tmp1->dimGlobal(),
                                             NULL);//tmp1->componentsNames());
      }
      else {
        result = new VectorSpace<V,M>(tmp2->env(),
                                             tmp2->prefix().c_str(),
                                             tmp2->dimGlobal(),
                                             NULL);//tmp2->componentsNames());
      }
    }
  }

  if (result == NULL) {
    const VectorSubset<V,M>* tmp1 = dynamic_cast<const VectorSubset<V,M>* >(&domain1);
    const VectorSubset<V,M>* tmp2 = dynamic_cast<const VectorSubset<V,M>* >(&domain2);

    if ((tmp1 != NULL) && (tmp2 != NULL)) {
      if (dim1 == dim2) {
        const BoxSubset<V,M>* box1 = dynamic_cast<const BoxSubset<V,M>* >(&domain1);
        const BoxSubset<V,M>* box2 = dynamic_cast<const BoxSubset<V,M>* >(&domain2);

        if ((box1 != NULL) && (box2 != NULL)) {
          V minV(box1->minValues());
          V maxV(box1->maxValues());
          for (unsigned int i = 0; i < dim1; ++i) {
            minV[i] = std::max(box1->minValues()[i],
                               box2->minValues()[i]);
          }
          for (unsigned int i = 0; i < dim1; ++i) {
            maxV[i] = std::min(box1->maxValues()[i],
                               box2->maxValues()[i]);
          }
          result = new BoxSubset<V,M>(box1->prefix().c_str(),
                                             box1->vectorSpace(),
                                             minV,
                                             maxV);
        }
        else {
          result = new IntersectionSubset<V,M>(tmp1->prefix().c_str(),
                                                      tmp1->vectorSpace(),
                                                      0., // FIX ME
                                                      domain1,
                                                      domain2);
        }
      }
      else {
        queso_error_msg("situation 001");
      }
    }
  }

  if (result == NULL) {
    const VectorSubset<V,M>* tmp1 = dynamic_cast<const VectorSubset<V,M>* >(&domain1);
    const VectorSpace <V,M>* tmp2 = dynamic_cast<const VectorSpace <V,M>* >(&domain2);

    if ((tmp1 != NULL) && (tmp2 != NULL)) {
      if (dim1 == dim2) {
        const BoxSubset<V,M>* box1 = dynamic_cast<const BoxSubset<V,M>* >(&domain1);
        if (box1 != NULL) {
          result = new BoxSubset<V,M>(box1->prefix().c_str(),
                                             box1->vectorSpace(),
                                             box1->minValues(),
                                             box1->maxValues());
        }
        else {
          queso_error_msg("situation 002");
        }
      }
      else {
        queso_error_msg("situation 003");
      }
    }
  }

  if (result == NULL) {
    const VectorSpace <V,M>* tmp1 = dynamic_cast<const VectorSpace <V,M>* >(&domain1);
    const VectorSubset<V,M>* tmp2 = dynamic_cast<const VectorSubset<V,M>* >(&domain2);

    if ((tmp1 != NULL) && (tmp2 != NULL)) {
      if (dim1 == dim2) {
        const BoxSubset<V,M>* box2 = dynamic_cast<const BoxSubset<V,M>* >(&domain2);
        if (box2 != NULL) {
          result = new BoxSubset<V,M>(box2->prefix().c_str(),
                                             box2->vectorSpace(),
                                             box2->minValues(),
                                             box2->maxValues());
        }
        else {
          queso_error_msg("situation 004");
        }
      }
      else {
        queso_error_msg("situation 005");
      }
    }
  }

  if (result == NULL) {
    queso_error_msg("situation 006");
  }

  return result;
}

}  // End namespace QUESO

template QUESO::VectorSet<QUESO::GslVector, QUESO::GslMatrix> *
QUESO::InstantiateIntersection(
    const QUESO::VectorSet<QUESO::GslVector, QUESO::GslMatrix>& domain1,
    const QUESO::VectorSet<QUESO::GslVector, QUESO::GslMatrix>& domain2);
