program copa__test_linear

  use copa__config, only : wp
  use copa__sampler, only : run_sampler
  use copa__store, only : store_chains
  use evortran__prng_rand, only : initialize_rands
  use evortran__prng_rand, only : randfloat
  use evortran__evolutions_float, only : evolve_population
  use evortran__individuals_float, only : individual

  implicit none

  integer :: i
  real(wp) :: data(20)
  real(wp) :: sigma_obs(20)
  real(wp) :: rand1
  real(wp) :: rand2
  real(wp), allocatable :: walkers(:,:)
  real(wp), allocatable :: chains(:,:,:)
  real(wp) :: ranges(2,3)
  type(individual) :: best_ind

  call initialize_rands(mode='twister', seed=1)

  ! Fake observed data for testing: y = 2*x + 1 + noise
  do i = 1, size(data)
    rand1 = randfloat(0.99e0_wp, 1.01e0_wp)
    rand2 = randfloat(0.99e0_wp, 1.01e0_wp)
    data(i) = rand1 * (2.0e0_wp * i + 10.0e0_wp - 4.0 * i ** 2) + rand2
  end do
  sigma_obs = 0.01e0_wp * data

  ranges(1,:) = -1.0e2_wp
  ranges(2,:) = 1.0e2_wp

  call run_sampler(  &
    3, log_prior, log_like,  &
    nsteps=100000,  &
    ranges=ranges,  &
    walkers=walkers,  &
    chains=chains)

  call store_chains(chains, "plots/linear/chains.npy", mode='machine')

  best_ind = evolve_population(  &
    1000, 3, fit_func,  &
    max_generations=1000, &
    selection_size=200,  &
    mating='sbx', &
    verbose=.true.,  &
    mutate='gaussian',  &
    lower_lim=-1.0e2_wp,  &
    upper_lim=1.0e2_wp)
  open(unit=11, file='plots/linear/evortran_best.csv', status='replace')
    do i = 1, size(best_ind%genes)
      write(11,'(F0.6)', advance='no') best_ind%genes(i)
      if (i < size(best_ind%genes)) then
        write(11,'(a)', advance='no') ','
      end if
    end do
  close(11)

contains

  subroutine log_prior(theta, logp)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(out) :: logp

    logp = -0.5_wp * sum((theta / 5.0_wp) ** 2)

  end subroutine log_prior

  subroutine log_like(theta, logl)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(out) :: logl

    real(wp) :: model_pred(size(data))
    integer :: x

    model_pred = [  &
      (theta(1) * real(x, wp) + theta(2) + theta(3) * real(x, wp) ** 2,  &
      x = 1, size(data))]
    logl = -0.5e0_wp * sum(((data - model_pred) / sigma_obs)**2)

  end subroutine log_like

  subroutine fit_func(ind, f)

    class(individual), intent(in) :: ind
    real(wp), intent(out) :: f

    call log_like(ind%genes, f)
    f = -f

  end subroutine fit_func

end program copa__test_linear
