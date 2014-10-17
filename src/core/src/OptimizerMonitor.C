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

#include <string>
#include <iomanip>
#include <queso/OptimizerMonitor.h>

namespace QUESO
{

  OptimizerMonitor::OptimizerMonitor( unsigned int n_iters )
    : m_display_conv(false),
      m_print_xmin(false)
  {
    m_minimizer_hist.reserve(n_iters);
    m_objective_hist.reserve(n_iters);
    m_norm_hist.reserve(n_iters);
  }

  OptimizerMonitor::~OptimizerMonitor()
  {}

  void OptimizerMonitor::set_display_output( bool enable_output, bool print_xmin )
  {
    m_display_conv = enable_output;
    m_print_xmin = print_xmin;
  }

  void OptimizerMonitor::append( std::vector<double>& x_min, double objective, double norm )
  {
    m_minimizer_hist.push_back(x_min);
    m_objective_hist.push_back(objective);
    m_norm_hist.push_back(norm);

    // Print out to screen if the user set the option
    if( m_display_conv )
      {
        // If we are appending the first entry, print a nice header
        if( m_minimizer_hist.size() == 1 )
          {
            this->print_header(std::cout);
          }

        // We're assuming here the size of the array is the current iteration
        this->print_iteration(m_norm_hist.size(),std::cout);
      }

  }

  void OptimizerMonitor::print_header( std::ostream& output ) const
  {
    unsigned int width = 35;

    output.width(5);
    output << "i";

    output.width(8);
    output << "f" << std::string(5,' ');

    output.width(11);
    output << "norm" << std::endl;
    output << std::string(width,'-') << std::endl;
  }

  void OptimizerMonitor::print_iteration( unsigned int iter, std::ostream& output ) const
  {
    output.width(5);
    output << iter;

    output.width(2);
    output << "  ";
    output.width(12);
    output << std::scientific << m_objective_hist[iter-1];

    output.width(2);
    output << "  ";
    output.width(12);
    output << m_norm_hist[iter-1] << std::endl;
  }

} // end namespace QUESO
