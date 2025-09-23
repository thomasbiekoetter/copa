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

  integer, parameter, public :: nwalkers_default = 20
  integer, parameter, public :: nsteps_default = 1000

end module copa__config
