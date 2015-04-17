!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! QUESO - a library to support the Quantification of Uncertainty
!! for Estimation, Simulation and Optimization
!!
!! Copyright (C) 2008,2009,2010 The PECOS Development Team
!!
!! This library is free software; you can redistribute it and/or
!! modify it under the terms of the Version 2.1 GNU Lesser General
!! Public License as published by the Free Software Foundation.
!!
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with this library; if not, write to the Free Software
!! Foundation, Inc. 51 Franklin Street, Fifth Floor,
!! Boston, MA  02110-1301  USA
!!
!!-----------------------------------------------------------------------el-
!!
!! Example MPI likelihood function ijobn F90.  This routine is called
!! from a C++ wrapper in exStatisticalInverseProblem1_likelihood.h and
!! is intended to match the pure C++ example provided in
!! examples/statisticalInverseProblem1.
!!
!!--------------------------------------------------------------------------

function f_likelihood(num_params,parameter_values,parameter_means,matrix) bind (C,name='f_likelihood')
  use iso_c_binding
  implicit none

  ! ------------------
  ! Incoming variables
  ! ------------------

  integer (C_int), value, intent(in) :: num_params                    ! total number of unknown parameters
  real (C_double),        intent(in) :: parameter_values (num_params) ! current parameter vector
  real (C_double),        intent(in) :: parameter_means  (num_params) ! mean parameter values
  real (C_double),        intent(in) :: matrix(num_params*num_params) ! matrix

  ! ------------------
  ! Return value
  ! ------------------

  real (C_double)      :: f_likelihood

  ! ------------------
  ! Local variables
  ! ------------------

  integer              :: irow,icol, istat
  integer, save        :: first_entry = 1
  real*8, allocatable  :: vector(:)

  ! ---------------------------------------------
  ! Initialization to perform on first entry only
  ! ---------------------------------------------

  if(first_entry .eq. 1)then
     allocate(vector(num_params),stat=istat)

     if(istat .ne. 0)then
        write(*,*) 'Error: Unable to allocate memory for vector'
        stop
     endif
  endif

  ! Compute likelihood (mimics C++ example in exStatisticalInverseProblem1_likelihood.h)

  f_likelihood = 0.0d0
  vector       = 0.0d0

  do irow = 1,num_params
     do icol = 1,num_params
        vector(irow) = vector(irow) + matrix(icol+num_params*(irow-1))*parameter_values(icol)
     enddo
  enddo

  do irow = 1,num_params
     f_likelihood = f_likelihood + parameter_values(irow)*vector(irow)
  enddo

end function f_likelihood
