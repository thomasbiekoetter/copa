module copa__parallel_sampler

  use copa__config, only : wp
  use copa__config, only : nwalkers_default
  use copa__config, only : nsteps_default
  use evortran__prng_rand, only : initialize_rands
  use evortran__prng_rand, only : randfloat
  use omp_lib

  implicit none

  private

  public :: run_parallel_sampler

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

  subroutine run_parallel_sampler(  &
    ndim, log_prior, log_like, nwalkers, nsteps, nthreads,  &
    ranges, walkers, chains)

    integer, intent(in) :: ndim
    procedure(log_prior_abstract) :: log_prior
    procedure(log_like_abstract) :: log_like
    integer, intent(in), optional :: nwalkers
    integer, intent(in), optional :: nsteps
    integer, intent(in), optional :: nthreads
    real(wp), intent(in), optional :: ranges(:,:)
    real(wp), intent(out), allocatable, optional :: walkers(:,:,:)
    real(wp), intent(out), allocatable, optional :: chains(:,:,:,:)

    integer :: nwal
    integer :: nste
    real(wp), allocatable :: wal(:,:,:)
    real(wp), allocatable :: cha(:,:,:,:)
    real(wp), allocatable :: ran(:,:)

    integer :: i
    integer :: j
    integer :: k
    integer :: step
    integer :: skip
    real(wp) :: z
    real(wp) :: q
    real(wp) :: log_p_current
    real(wp) :: log_p_proposed
    real(wp) :: rand
    real(wp) :: new_pos(ndim)
    real(wp), parameter :: a = 2.0e0_wp
    integer :: ncpu
    integer :: nthr

    if (present(nwalkers)) then
      nwal = nwalkers
    else
      nwal = nwalkers_default
    end if

    if (present(nsteps)) then
      nste = nsteps
    else
      nste = nsteps_default
    end if

    ncpu = omp_get_num_procs()
    if (present(nthreads)) then
      nthr = nthreads
    else
      nthr = ncpu
    end if
    call omp_set_num_threads(nthr)
    write(*,'(A,I0,A,I0,A)')  &
      "Running with ", nthr, " threads out of ", ncpu, " available CPUs."

    if (present(ranges)) then
      if (.not. (rank(ranges) == 2)) then
        write(*,*) "Optional argument 'ranges' should have rank 2."
        call exit
      end if
      if (.not. (size(ranges, 1) == 2) .and. (size(ranges, 2) == ndim)) then
        write(*,*) "Optional argument 'ranges' should have shape (2, ndim)."
        call exit
      end if
      ran = ranges
    else
      allocate(ran(2,ndim))
      ran(1,:) = 0.0e0_wp
      ran(2,:) = 1.0e0_wp
    end if

    allocate(wal(ndim, nwal, nthr))
    allocate(cha(ndim, nwal, nste, nthr))

    !$omp parallel do  &
    !$omp default(none)  &
    !$omp private(  &
    !$omp   k, i, j, step, rand, z, new_pos,  &
    !$omp   log_p_current, log_p_proposed, q, skip)  &
    !$omp shared(  &
    !$omp   nthr, ndim, nwal, wal, ran, nste, cha)
    do k = 1, nthr

      do i = 1, ndim
        do j = 1, nwal
          wal(i,j,k) = randfloat(ran(1,i), ran(2,i))
        end do
      end do

      do step = 1, nste

        do i = 1, nwal

          j = i
          do while(j == i)
            rand = randfloat()
            j = int(rand * nwal) + 1
          end do

          ! Stretch factor z ~ 1/sqrt(z)
          rand = randfloat()
          z = ((a - 1.0e0_wp / a) * rand + 1.0e0_wp / a)
          z = z**2

          ! Propose new position
          new_pos = wal(:,j,k) + z * (wal(:,i,k) - wal(:,j,k))

          ! Log probabilities
          call log_prob(wal(:,i,k), log_p_current)
          call log_prob(new_pos, log_p_proposed)

          q = z**(ndim - 1) * exp(log_p_proposed - log_p_current)

          rand = randfloat()
          if (rand < min(1.0e0_wp, q)) then
            wal(:,i,k) = new_pos
          end if

        end do

        ! Store current state
        cha(:,:,step,k) = wal(:,:,k)

        ! Print progress
        if (k == 1) then
          skip = int(nste  / 10)
          if (mod(step, skip) == 0) then
            write(*,'(a,i8,a,*(f8.4,1x))') 'Step:', step, '  Walker 1:', wal(:,1,k)
          end if
        end if

      end do

    end do
    !$omp end parallel do

    if (present(walkers)) then
      walkers = wal
    end if

    if (present(chains)) then
      chains = cha
    end if

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

  end subroutine run_parallel_sampler

end module copa__parallel_sampler
