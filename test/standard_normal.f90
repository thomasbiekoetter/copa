program copa__test_standard_normal

  use copa__config, only : wp
  use copa__sampler, only : run_sampler
  use copa__store, only : store_chains
  use copa__prior_functions, only : uniform_prior
  use evortran__prng_rand, only : initialize_rands

  implicit none

  integer, parameter :: ndim = 3
  real(wp), parameter :: lower(ndim) = -10.0e0_wp
  real(wp), parameter :: upper(ndim) = 10.0e0_wp

  real(wp), allocatable :: walkers(:,:)
  real(wp), allocatable :: chains(:,:,:)
  real(wp) :: ranges(2,ndim)

  ranges(1,:) = lower(1)
  ranges(2,:) = upper(1)

  call initialize_rands(mode='twister', seed=1)

  call run_sampler(  &
    ndim, log_prior, log_like,  &
    nsteps=10000,  &
    ranges=ranges,  &
    walkers=walkers,  &
    chains=chains)

  call store_chains(  &
    chains,  &
    "plots/standard_normal/chains.npy",  &
    mode='machine')

contains

  subroutine log_prior(theta, logp)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(out) :: logp

    call uniform_prior(theta, lower, upper, logp)

  end subroutine log_prior

  subroutine log_like(theta, logl)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(out) :: logl

    logl = -0.5e0_wp * sum(theta ** 2)

  end subroutine log_like

end program copa__test_standard_normal
