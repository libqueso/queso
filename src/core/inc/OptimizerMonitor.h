//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
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

#ifndef UQ_OPTIMIZER_MONITOR_H
#define UQ_OPTIMIZER_MONITOR_H

#include <vector>
#include <iostream>

namespace QUESO
{
  //! Object to monitor convergence of optimizers
  class OptimizerMonitor
  {
  public:

    //! Constructor
    /*!
     * Input paramter n_iters will reserve space for monitor data
     * for n_iters. This is for efficiency - if the number of iterations
     * goes beyond n_iters, there will not be an error, just less performance
     * by the monitor. Defaults to reserving space for 100 iterations
     */
    OptimizerMonitor( unsigned int n_iters = 100 );

    ~OptimizerMonitor();

    void set_display_output( bool enable_output, bool print_xmin );
    
    void append( std::vector<double>& x_min, double objective, double norm );
    
    void reset();

  private:
    
    bool m_display_conv;
    bool m_print_xmin;

    std::vector<std::vector<double> > m_minimizer_hist;
    std::vector<double> m_objective_hist;
    std::vector<double> m_norm_hist;

    void print_header( std::ostream& output ) const;
    void print_iteration( unsigned int iter, std::ostream& output ) const;

  };

} // end namespace QUESO

#endif // UQ_OPTIMIZER_MONITOR_H
