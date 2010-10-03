//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef __UQ_INSTANTIATE_INTERSECTION_H__
#define __UQ_INSTANTIATE_INTERSECTION_H__

#include <uqVectorSpace.h>
#include <uqVectorSubset.h>

template<class V, class M>
uqVectorSetClass<V,M>*
uqInstantiateIntersection(const uqVectorSetClass<V,M>& domain1, const uqVectorSetClass<V,M>& domain2)
{
  uqVectorSetClass<V,M>* result = NULL;

  unsigned int dim1 = domain1.vectorSpace().dimGlobal();
  unsigned int dim2 = domain2.vectorSpace().dimGlobal();

  if (result == NULL) {
    const uqVectorSpaceClass<V,M>* tmp1 = dynamic_cast<const uqVectorSpaceClass<V,M>* >(&domain1);
    const uqVectorSpaceClass<V,M>* tmp2 = dynamic_cast<const uqVectorSpaceClass<V,M>* >(&domain2);

    if ((tmp1 != NULL) && (tmp2 != NULL)) {
      if (dim1 < dim2) {
        result = new uqVectorSpaceClass<V,M>(tmp1->env(),
                                             tmp1->prefix().c_str(),
                                             tmp1->dimGlobal(),
                                             NULL);//tmp1->componentsNames());
      }
      else if (dim1 == dim2) {
        result = new uqVectorSpaceClass<V,M>(tmp1->env(),
                                             tmp1->prefix().c_str(),
                                             tmp1->dimGlobal(),
                                             NULL);//tmp1->componentsNames());
      }
      else {
        result = new uqVectorSpaceClass<V,M>(tmp2->env(),
                                             tmp2->prefix().c_str(),
                                             tmp2->dimGlobal(),
                                             NULL);//tmp2->componentsNames());
      }
    }
  }

  if (result == NULL) {
    const uqVectorSubsetClass<V,M>* tmp1 = dynamic_cast<const uqVectorSubsetClass<V,M>* >(&domain1);
    const uqVectorSubsetClass<V,M>* tmp2 = dynamic_cast<const uqVectorSubsetClass<V,M>* >(&domain2);

    if ((tmp1 != NULL) && (tmp2 != NULL)) {
      if (dim1 == dim2) {
        const uqBoxSubsetClass<V,M>* box1 = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&domain1);
        const uqBoxSubsetClass<V,M>* box2 = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&domain2);

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
          result = new uqBoxSubsetClass<V,M>(box1->prefix().c_str(),
                                             box1->vectorSpace(),
                                             minV,
                                             maxV);
        }
        else {
          result = new uqIntersectionSubsetClass<V,M>(tmp1->prefix().c_str(),
                                                      tmp1->vectorSpace(),
                                                      0., // FIX ME
                                                      domain1,
                                                      domain2);
        }
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            domain1.env().worldRank(),
                            "uqInstantiateIntersection<V,M>()",
                            "situation 001");
      }
    }
  }

  if (result == NULL) {
    const uqVectorSubsetClass<V,M>* tmp1 = dynamic_cast<const uqVectorSubsetClass<V,M>* >(&domain1);
    const uqVectorSpaceClass <V,M>* tmp2 = dynamic_cast<const uqVectorSpaceClass <V,M>* >(&domain2);

    if ((tmp1 != NULL) && (tmp2 != NULL)) {
      if (dim1 == dim2) {
        const uqBoxSubsetClass<V,M>* box1 = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&domain1);
        if (box1 != NULL) {
          result = new uqBoxSubsetClass<V,M>(box1->prefix().c_str(),
                                             box1->vectorSpace(),
                                             box1->minValues(),
                                             box1->maxValues());
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              domain1.env().worldRank(),
                              "uqInstantiateIntersection<V,M>()",
                              "situation 002");
        }
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            domain1.env().worldRank(),
                            "uqInstantiateIntersection<V,M>()",
                            "situation 003");
      }
    }
  }

  if (result == NULL) {
    const uqVectorSpaceClass <V,M>* tmp1 = dynamic_cast<const uqVectorSpaceClass <V,M>* >(&domain1);
    const uqVectorSubsetClass<V,M>* tmp2 = dynamic_cast<const uqVectorSubsetClass<V,M>* >(&domain2);

    if ((tmp1 != NULL) && (tmp2 != NULL)) {
      if (dim1 == dim2) {
        const uqBoxSubsetClass<V,M>* box2 = dynamic_cast<const uqBoxSubsetClass<V,M>* >(&domain2);
        if (box2 != NULL) {
          result = new uqBoxSubsetClass<V,M>(box2->prefix().c_str(),
                                             box2->vectorSpace(),
                                             box2->minValues(),
                                             box2->maxValues());
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              domain1.env().worldRank(),
                              "uqInstantiateIntersection<V,M>()",
                              "situation 004");
        }
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            domain1.env().worldRank(),
                            "uqInstantiateIntersection<V,M>()",
                            "situation 005");
      }
    }
  }

  if (result == NULL) {
    UQ_FATAL_TEST_MACRO(true,
                        domain1.env().worldRank(),
                        "uqInstantiateIntersection<V,M>()",
                        "situation 006");
  }

  return result;
}
#endif // __UQ_INSTANTIATE_INTERSECTION_H__

