module copa__config

  implicit none

  private

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: qp = selected_real_kind(30,4931)

#ifdef QUAD
  integer, parameter, public :: wp = qp
#else
  integer, parameter, public :: wp = dp
#endif

  integer, parameter, public :: nwalkers_default = 100
  integer, parameter, public :: nsteps_default = 1000
  integer, parameter, public :: nensembles_default = 4
  real(wp), parameter, public :: a_default = 2.0e0_wp
  character(len=*), parameter, public :: parallel_method_default = 'redblack'

end module copa__config
