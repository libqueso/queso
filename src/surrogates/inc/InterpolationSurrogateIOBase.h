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

#ifndef UQ_INTERPOLATION_SURROGATE_IO_BASE_H
#define UQ_INTERPOLATION_SURROGATE_IO_BASE_H

// QUESO
#include <queso/ScopedPtr.h>
#include <queso/InterpolationSurrogateData.h>

namespace QUESO
{
  template<class V, class M>
  class InterpolationSurrogateIOBase
  {
  public:

    InterpolationSurrogateIOBase();

    virtual ~InterpolationSurrogateIOBase(){};

    virtual void read( const std::string& filename,
                       const FullEnvironment& env,
                       const std::string& vector_space_prefix,
                       int reading_rank ) =0;

    virtual void write( const std::string& filename,
                        const InterpolationSurrogateData<V,M>& data,
                        int writing_rank ) const =0;

    //! Reference to data object
    /*! If we read the data, then use this method to grab
        a reference to the data in order to construct
        an InterpolationSurrogateBase subclass. */
    const InterpolationSurrogateData<V,M>& data() const
    { return *(this->m_data.get()); };

  protected:

    typename ScopedPtr<VectorSpace<V,M> >::Type m_vector_space;

    typename ScopedPtr<BoxSubset<V,M> >::Type m_domain;

    std::vector<unsigned int> m_n_points;

    typename ScopedPtr<InterpolationSurrogateData<V,M> >::Type m_data;
  };

} // end namespace QUESO

#endif // UQ_INTERPOLATION_SURROGATE_IO_BASE_H
