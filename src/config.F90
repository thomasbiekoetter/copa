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
  character(len=*), parameter, public :: store_mode_default = 'machine'
  logical, parameter, public :: store_separate_default = .true.

end module copa__config
