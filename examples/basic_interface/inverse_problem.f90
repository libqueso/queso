
!----------------------------------------
! A simple Calibration example using the
! basic QUESO interface..
!---------------------------------------

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

