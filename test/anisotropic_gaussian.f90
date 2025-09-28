program copa__test_anisotropic_gaussian

  use copa__config, only : wp
  use copa__sampler, only : run_sampler
  use copa__store, only : store_chains
  use copa__prior_functions, only : uniform_prior
  use evortran__prng_rand, only : initialize_rands

  implicit none

  integer, parameter :: ndim = 5
  integer, parameter :: nsteps = 10000
  real(wp), parameter :: lower(ndim) = -1.0e2_wp
  real(wp), parameter :: upper(ndim) = 1.0e2_wp
  real(wp) :: mu(ndim)
  real(wp) :: var(ndim)
  real(wp) :: sigma_inv(ndim, ndim)

  real(wp), allocatable :: walkers(:,:)
  real(wp), allocatable :: chains(:,:,:)
  real(wp) :: ranges(2,ndim)

  integer :: i

  sigma_inv = 0.0e0_wp
  do i = 1, ndim
    mu(i) = real(i, wp)
    var(i) = mu(i) * 1.0e-1_wp
    sigma_inv(i,i) = 1.0e0_wp / var(i)
  end do

  ranges(1,:) = lower(1)
  ranges(2,:) = upper(1)

  call initialize_rands(mode='twister', seed=1)

  call run_sampler(  &
    ndim, log_prior, log_like,  &
    nsteps=nsteps,  &
    ranges=ranges,  &
    walkers=walkers,  &
    chains=chains)

  call store_chains(  &
    chains,  &
    "plots/anisotropic_gaussian/chains.npy",  &
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

    real(wp) :: x(size(theta))

    x = theta - mu
    logl = -0.5e0_wp * dot_product(x, matmul(sigma_inv, x))

  end subroutine log_like

end program copa__test_anisotropic_gaussian
