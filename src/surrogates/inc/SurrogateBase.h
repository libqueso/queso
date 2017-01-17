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

#ifndef UQ_SURROGATE_BASE_H
#define UQ_SURROGATE_BASE_H

namespace QUESO
{
  class GslVector;

  //! Base class for surrogates of models
  /*! Defines basic interface for using surrogates of models. These surrogates
      map an \f$ n\f$ dimensional parameter space to the reals. That is
      \f$ f: \mathbb{R}^n \rightarrow \mathbb{R} \f$. The idea
      is that we have some surrogate of the parameter-to-data map so that
      it can be used in the likelihood classes. Subclasses will define the
      particular surrogate model. Other classes will be used to build up the
      surrogate from the user's model. */
  template<class V = GslVector>
  class SurrogateBase
  {
  public:

    SurrogateBase(){};

    virtual ~SurrogateBase(){};

    //! Method to return value given the parameter vector
    virtual double evaluate(const V & domainVector) const =0;

  };

} // end namespace QUESO

#endif // UQ_SURROGATE_BASE_H
