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
!! Example likelihood function in F90.  This routine is called from a
!! C++ wrapper in exStatisticalInverseProblem1_likelihood.h and is
!! intended to match the pure C++ example provided in
!! examples/statisticalInverseProblem1.
!!
!!--------------------------------------------------------------------------

function f_likelihood(num_params,parameter_values,parameter_means,matrix,subid,subcomm) bind (C,name='f_likelihood')
  use iso_c_binding
  implicit none
  include "mpif.h"

  ! ------------------
  ! Incoming variables
  ! ------------------

  integer (C_int), value, intent(in) :: num_params                    ! total number of unknown parameters
  real    (C_double),     intent(in) :: parameter_values (num_params) ! current parameter vector
  real    (C_double),     intent(in) :: parameter_means  (num_params) ! mean parameter values
  real    (C_double),     intent(in) :: matrix(num_params*num_params) ! matrix
  integer (C_int),value,  intent(in) :: subid                         ! QUESO subcommunicator ID
  integer (C_int),value,  intent(in) :: subcomm                       ! MPI subcommunicator for doing work

  ! ------------------
  ! Return value
  ! ------------------

  real (C_double)      :: f_likelihood

  ! ------------------
  ! Local variables
  ! ------------------

  integer                   :: irow,icol, istat, ierr
  integer                   :: num_global_procs, num_subcomm_procs
  integer                   :: rank_global, rank_subcomm
  integer, save             :: first_entry = 1
  real*8, allocatable, save :: vector(:)

  integer, save :: foo = 0

  ! ---------------------------------------------
  ! Initialization to perform on first entry only
  ! ---------------------------------------------

  if(first_entry .eq. 1)then
     allocate(vector(num_params),stat=istat)

     if(istat .ne. 0)then
        write(*,*) 'Error: Unable to allocate memory for vector'
        stop
     endif

     ! determine global/local MPI environment info

     call MPI_Comm_size(MPI_COMM_WORLD,num_global_procs,ierr)
     call MPI_Comm_rank(MPI_COMM_WORLD,rank_global,ierr)

     call MPI_Comm_size(subcomm,num_subcomm_procs,ierr)
     call MPI_Comm_rank(subcomm,rank_subcomm,ierr)

     if(rank_global .eq. 0)then
        write(*,*) ' '
        write(*,*) 'Likelihood parallel environment (F90):'
        write(*,*) '--> Total # of MPI processors available     = ', num_global_procs
        write(*,*) '--> Total # of SubEnvironments requested    = ', num_global_procs/num_subcomm_procs
        write(*,*) '--> Number of MPI tasks per SubEnvironment  = ', num_subcomm_procs
        write(*,*) ' '
     endif

     first_entry = 0
  endif

  ! -------------------------------------------------------------------
  ! Compute likelihood (mimics C++ example in
  ! exStatisticalInverseProblem1_likelihood.h).
  !
  ! In practical applications, this would be the place to perform a
  ! forward problem evaluation to determine the likelihood.
  ! -------------------------------------------------------------------

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

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

end function f_likelihood
