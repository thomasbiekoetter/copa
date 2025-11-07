program copa__test_rosenbrock

  use copa__config, only : wp
  use copa__parallel_sampler, only : run_parallel_sampler
  use copa__store, only : store_chains
  use copa__store, only : store_log_probs
  use copa__prior_functions, only : uniform_prior
  use evortran__prng_rand, only : initialize_rands

  implicit none

  real(wp), parameter :: a = 1.0e0_wp
  real(wp), parameter :: b = 1.0e2_wp

  integer, parameter :: ndim = 2
  integer, parameter :: nthreads = 8
  integer, parameter :: nsteps = 10000
  real(wp), parameter :: lower(ndim) = [  &
    -2.0e0_wp,  &
    -1.0e0_wp]
  real(wp), parameter :: upper(ndim) = [  &
    2.0e0_wp,  &
    3.0e0_wp]
  real(wp), allocatable :: walkers(:,:,:)
  real(wp), allocatable :: chains(:,:,:,:)
  real(wp), allocatable :: log_probs(:,:,:)
  real(wp) :: ranges(2,ndim)

  call initialize_rands(mode="twister", seed=1)

  ranges(1,:) = lower
  ranges(2,:) = upper

  call run_parallel_sampler(  &
    ndim, log_prior, log_like,  &
    nsteps=nsteps,  &
    nthreads=nthreads,  &
    ranges=ranges,  &
    walkers=walkers,  &
    chains=chains,  &
    log_probs=log_probs)

  call store_chains(  &
    chains,  &
    "plots/rosenbrock/chains.npy",  &
    mode="machine")

  call store_log_probs(  &
    log_probs,  &
    "plots/rosenbrock/log_probs.npy",  &
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

    real(wp) :: x
    real(wp) :: y

    x = theta(1)
    y = theta(2)

    logl = -((a - x) ** 2 + b * (y - x ** 2) ** 2)

  end subroutine log_like

end program copa__test_rosenbrock
