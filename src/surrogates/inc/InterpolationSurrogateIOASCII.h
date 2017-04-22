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

#ifndef UQ_INTERPOLATION_SURROGATE_IO_ASCII_H
#define UQ_INTERPOLATION_SURROGATE_IO_ASCII_H

#include <queso/InterpolationSurrogateIOBase.h>

namespace QUESO
{
  template<class V, class M>
  class InterpolationSurrogateIOASCII : public InterpolationSurrogateIOBase<V,M>
  {
  public:

    InterpolationSurrogateIOASCII();

    virtual ~InterpolationSurrogateIOASCII(){};

    //! Read Interpolation surrogate data from filename using processor reading_rank
    /*! This will read the data the file given by filename and setup all
        the infrastructure for InterpolationSurrogateData so that the user
        can then call the data() method. This can then be used to construct
        an Interpolation object. env.fullRank() must contain reading_rank.
        By default, processor 0 reads the data. */
    virtual void read( const std::string& filename,
                       const FullEnvironment& env,
                       const std::string& vector_space_prefix,
                       int reading_rank = 0 );

    //! Write interpolation surrogate data to filename using processor writing_rank
    /*! env.fullRank() must contain writing_rank. By default processor 0
        writes the data. */
    virtual void write( const std::string& filename,
                        const InterpolationSurrogateData<V,M>& data,
                        int writing_rank = 0 ) const;

  };
} // end namespace QUESO

#endif // UQ_INTERPOLATION_SURROGATE_IO_ASCII_H
