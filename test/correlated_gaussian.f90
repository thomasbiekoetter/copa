program copa__test_correlated_gaussian

  use copa__config, only : wp
  use copa__sampler, only : run_sampler
  use copa__store, only : store_chains
  use copa__prior_functions, only : uniform_prior
  use evortran__prng_rand, only : initialize_rands
  use linalg_matrices_inverse, only : inversereal
  use linalg_matrices_diagonalize, only : diagonalizerealsymmetric

  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nsteps = 10000
  real(wp), parameter :: lower(ndim) = -1.0e2_wp
  real(wp), parameter :: upper(ndim) = 1.0e2_wp
  real(wp) :: mu(ndim)
  real(wp) :: sigma(ndim,ndim)
  real(wp) :: sigma_inv(ndim,ndim)

  real(wp), allocatable :: walkers(:,:)
  real(wp), allocatable :: chains(:,:,:)
  real(wp) :: ranges(2,ndim)

  mu = 0.0e0_wp

  sigma = reshape([  &
    1.0e0_wp,  0.8e0_wp,  0.1e0_wp,  &
    0.8e0_wp,  1.0e0_wp, -0.4e0_wp,  &
    0.1e0_wp, -0.4e0_wp,  1.0e0_wp],  &
    [ndim, ndim])

  call check_model()

  sigma_inv = inversereal(sigma)

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
    "plots/correlated_gaussian/chains.npy",  &
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

  subroutine check_model()

    real(wp) :: sigma_eigval(ndim)
    real(wp) :: sigma_eigvec(ndim,ndim)
    integer :: order

    order = 1
    call diagonalizerealsymmetric(  &
      sigma, sigma_eigval, sigma_eigvec, order)

    if (abs(sum(transpose(sigma) - sigma)) > 1.0e-10_wp) then
      write(*,*) "Error: correlation matrix not symmetric."
      call exit
    end if

    if (any(sigma_eigval < 0.0e0_wp)) then
      write(*,*) "Error: correlation matrix not positive definite."
      call exit
    end if

  end subroutine check_model

end program copa__test_correlated_gaussian
