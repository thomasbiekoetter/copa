module copa__redblack_parallel_sampler

  use copa__config, only : wp
  use copa__config, only : nensembles_default
  use copa__config, only : nwalkers_default
  use copa__config, only : nsteps_default
  use copa__config, only : a_default
  use evortran__prng_rand, only : are_rands_initialized
  use evortran__prng_rand, only : initialize_rands
  use evortran__prng_rand, only : randfloat
  use iso_fortran_env, only : error_unit
  use omp_lib

  implicit none

  private

  public :: run_redblack_parallel_sampler

  abstract interface

    subroutine log_prior_abstract(theta, logp)
      import :: wp
      implicit none
      real(wp), intent(in) :: theta(:)
      real(wp), intent(out) :: logp
    end subroutine log_prior_abstract

    subroutine log_like_abstract(theta, logl)
      import :: wp
      implicit none
      real(wp), intent(in) :: theta(:)
      real(wp), intent(out) :: logl
    end subroutine log_like_abstract

  end interface

contains

  subroutine run_redblack_parallel_sampler(  &
    ndim, log_prior, log_like,  &
    nwalkers, nsteps, nthreads, nensembles,  &
    ranges, walkers, chains, log_probs)
    integer, intent(in) :: ndim
    procedure(log_prior_abstract) :: log_prior
    procedure(log_like_abstract) :: log_like
    integer, intent(in), optional :: nwalkers
    integer, intent(in), optional :: nsteps
    integer, intent(in), optional :: nthreads
    integer, intent(in), optional :: nensembles
    real(wp), intent(in), optional :: ranges(:,:)
    real(wp), intent(out), allocatable, optional :: walkers(:,:,:)
    real(wp), intent(out), allocatable, optional :: chains(:,:,:,:)
    real(wp), intent(out), allocatable, optional :: log_probs(:,:,:)

    integer :: nwal
    integer :: nste
    integer :: nens
    real(wp), allocatable :: wal(:, :, :)
    real(wp), allocatable :: cha(:, :, :, :)
    real(wp), allocatable :: ran(:, :)
    real(wp), allocatable :: lg_pb(:, :, :)
    real(wp), allocatable :: lp(:, :)

    integer :: i
    integer :: j
    integer :: e
    integer :: step
    real(wp) :: z
    real(wp) :: log_q
    real(wp) :: log_p_proposed
    real(wp) :: rand
    real(wp) :: new_pos(ndim)
    real(wp) :: a
    integer :: ncpu
    integer :: nthr
    integer :: nhalf
    integer :: half
    integer :: ia
    integer :: ib
    integer :: ja
    integer :: jb
    integer, allocatable :: naccept(:)
    integer :: naccept_half
    integer :: skip

    a = a_default

    nwal = nwalkers_default
    if (present(nwalkers)) nwal = nwalkers
    if (nwal <= 1) error stop "nwalkers must be >= 2."
    if (mod(nwal, 2) /= 0) error stop "nwalkers must be even."

    nste = nsteps_default
    if (present(nsteps)) nste = nsteps

    nens = nensembles_default
    if (present(nensembles)) nens = nensembles
    if (nens < 1) error stop "nensembles must be >= 1."

    ncpu = omp_get_num_procs()
    nthr = ncpu
    if (present(nthreads)) nthr = nthreads
    if (nthr < 2) error stop "nthreads must be >= 2."
    if (nthr > ncpu) &
      write(error_unit,'(a,i0,a,i0,a)') &
        "Warning: nthreads (", nthr, ") exceeds available CPUs (", ncpu, ")."
    call omp_set_num_threads(nthr)
    write(*,'(A,I0,A,I0,A)')  &
      "Running with ", nthr, " threads out of ", ncpu, " available CPUs."

    if (present(ranges)) then
      if (rank(ranges) /= 2) error stop  &
        "Optional argument 'ranges' should have rank 2."
      if ((size(ranges, 1) /= 2) .or. (size(ranges, 2) /= ndim)) error stop  &
        "Optional argument 'ranges' should have shape (2, ndim)."
      ran = ranges
    else
      allocate(ran(2, ndim))
      ran(1, :) = 0.0e0_wp
      ran(2, :) = 1.0e0_wp
    end if

    if (.not. are_rands_initialized)  &
      call initialize_rands(  &
        mode='twister', seed=1, nthreads=nthr)

    allocate(wal(nens, ndim, nwal))
    allocate(cha(nens, ndim, nwal, nste))
    allocate(lg_pb(nens, nwal, nste))
    allocate(lp(nens, nwal))

    allocate(naccept(nens))
    nhalf = nwal / 2
    naccept = 0
    skip = max(1, int(nste / 10))

    do e = 1, nens

      do j = 1, nwal
        do i = 1, ndim
          wal(e, i, j) = randfloat(ran(1, i), ran(2, i))
        end do
        call log_prob(wal(e, :, j), lp(e, j))
      end do

      do step = 1, nste

        do half = 0, 1

          if (half == 0) then
            ia = 1
            ib = nhalf
            ja = nhalf + 1
            jb = nwal
          else
            ia = nhalf + 1
            ib = nwal
            ja = 1
            jb = nhalf
          end if

          naccept_half = 0

          !$omp parallel do default(none) &
          !$omp private(i, j, rand, z, new_pos, log_p_proposed, log_q) &
          !$omp shared(a, e, ia, ib, ja, jb, ndim, wal, lp)  &
          !$omp reduction(+:naccept_half)
          do i = ia, ib

            ! partner drawn only from the frozen half
            rand = randfloat()
            j = ja + int(rand * (jb - ja + 1))
            j = min(j, jb)

            ! Stretch factor z ~ 1/sqrt(z)
            rand = randfloat()
            z = ((a - 1.0e0_wp) * rand + 1.0e0_wp)**2 / a

            ! Propose new position
            new_pos = wal(e, :, j) + z * (wal(e, :, i) - wal(e, :, j))

            ! Log probabilities
            call log_prob(new_pos, log_p_proposed)

!           q = z ** (ndim - 1) * exp(log_p_proposed - lp(e, i))
            log_q = real(ndim - 1, wp) * log(z) + (log_p_proposed - lp(e, i))

            rand = randfloat()
            if (log(rand) < log_q) then ! equivalent to rand < min(1, exp(log_q))
              wal(e, :, i) = new_pos
              lp(e, i) = log_p_proposed
              naccept_half = naccept_half + 1
            end if

          end do
          !$omp end parallel do

          naccept(e) = naccept(e) + naccept_half

        end do

        cha(e, :, :, step) = wal(e, :, :)
        lg_pb(e, :, step)  = lp(e, :)

        if (mod(step, skip) == 0) then
          write(*,'(a,i0,a,i0,a,i8,a,f8.3)') &
            'Ensemble ', e, '/', nens, '   Step:', step, &
            '   accept:', real(naccept(e), wp) / real(nwal * step * 2, wp)
        end if

      end do

    end do

    if (present(walkers)) walkers = wal
    if (present(chains)) chains = cha
    if (present(log_probs)) log_probs = lg_pb

    contains

      subroutine log_prob(theta, logpb)

        real(wp), intent(in) :: theta(:)
        real(wp), intent(out) :: logpb

        real(wp) :: logp
        real(wp) :: logl

        call log_prior(theta, logp)
        call log_like(theta, logl)

        logpb = logp + logl

      end subroutine log_prob

  end subroutine run_redblack_parallel_sampler

end module copa__redblack_parallel_sampler
