!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
! Copyright (C) 2008,2009 The PECOS Development Team
!
! Please see http://pecos.ices.utexas.edu for more information.
!
! This file is part of the QUESO Library (Quantification of Uncertainty
! for Estimation, Simulation and Optimization).
!
! QUESO is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! QUESO is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with QUESO. If not, see <http://www.gnu.org/licenses/>.
!
!--------------------------------------------------------------------------
!
! $Id$
!
! Simple toy problem to illustrate the use of the Basic QUESO
! interface to perform a statistical inverse problem in Fortran..
! 
!--------------------------------------------------------------------------
!-------------------------------------------------------------------------- 

program main
  integer :: ierr
  external my_likelihood

  call MPI_Init(ierr);

  print *,'--> Initializing QUESO Environment...'

  call QUESO_init('koomie.inp');                  ! initialize QUESO environment and define input file
  call QUESO_statistical_inversion(my_likelihood) ! register my likelihood function with QUESO
  call QUESO_finalize()                           ! Run MCMC  

  call MPI_Finalize(ierr);
  stop
end program main

real*8 function my_likelihood(params)
  implicit none
  
  real*8,intent(in)  :: params

  real*8 :: needle = 42.
  real*8 :: sigma  = 1.

  my_likelihood = (params - needle)**2/(sigma**2)

  return
end function my_likelihood

