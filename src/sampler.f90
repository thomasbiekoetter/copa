module copa__sampler

  use copa__config, only : wp
  use copa__config, only : nwalkers_default
  use copa__config, only : nsteps_default
  use evortran__prng_rand, only : initialize_rands
  use evortran__prng_rand, only : randfloat

  implicit none

  private

  public :: run_sampler

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

  subroutine run_sampler(  &
    ndim, log_prior, log_like,  &
    nwalkers, nsteps, ranges, walkers, chains, log_probs)

    integer, intent(in) :: ndim
    procedure(log_prior_abstract) :: log_prior
    procedure(log_like_abstract) :: log_like
    integer, intent(in), optional :: nwalkers
    integer, intent(in), optional :: nsteps
    real(wp), intent(in), optional :: ranges(:,:)
    real(wp), intent(out), allocatable, optional :: walkers(:,:)
    real(wp), intent(out), allocatable, optional :: chains(:,:,:)
    real(wp), intent(out), allocatable, optional :: log_probs(:,:)

    real(wp) :: x
    integer :: nwal
    integer :: nste
    real(wp), allocatable :: wal(:,:)
    real(wp), allocatable :: cha(:,:,:)
    real(wp), allocatable :: ran(:,:)
    real(wp), allocatable :: lg_pb(:,:)

    integer :: i
    integer :: j
    integer :: step
    integer :: skip
    real(wp) :: z
    real(wp) :: q
    real(wp) :: log_p_current
    real(wp) :: log_p_proposed
    real(wp) :: rand
    real(wp) :: new_pos(ndim)
    real(wp), parameter :: a = 2.0e0_wp

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

    allocate(wal(ndim, nwal))
    allocate(cha(ndim, nwal, nste))
    allocate(lg_pb(nwal,nste))

    do i = 1, ndim
      do j = 1, nwal
        wal(i,j) = randfloat(ran(1,i), ran(2,i))
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
        new_pos = wal(:,j) + z * (wal(:,i) - wal(:,j))

        ! Log probabilities
        call log_prob(wal(:,i), log_p_current)
        call log_prob(new_pos, log_p_proposed)

        q = z**(ndim - 1) * exp(log_p_proposed - log_p_current)

        rand = randfloat()
        if (rand < min(1.0e0_wp, q)) then
          wal(:,i) = new_pos
          lg_pb(i,step) = log_p_proposed
        else
          lg_pb(i,step) = log_p_current
        end if

      end do

      ! Store current state
      cha(:,:,step) = wal

      ! Print progress
      skip = int(nste  / 10)
      if (mod(step, skip) == 0) then
        write(*,'(a,i8,a,*(f8.4,1x))') 'Step:', step, '  Walker 1:', wal(:,1)
      end if

    end do

    if (present(walkers)) then
      walkers = wal
    end if

    if (present(chains)) then
      chains = cha
    end if

    if (present(log_probs)) then
      log_probs = lg_pb
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

  end subroutine run_sampler

end module copa__sampler
