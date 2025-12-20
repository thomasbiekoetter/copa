module copa__store

  use copa__config, only : wp

  implicit none

  private :: store_chains_human_single
  private :: store_chains_human_parallel
  private :: store_chains_machine_single
  private :: store_chains_machine_parallel

  private :: store_log_probs_human_single
  private :: store_log_probs_human_parallel
  private :: store_log_probs_machine_single
  private :: store_log_probs_machine_parallel

  public :: store_chains
  public :: store_log_probs

contains

  subroutine store_chains(chains, filename, mode)

    real(wp), intent(in) :: chains(..)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: mode

    character(len=:), allocatable :: md

    if (present(mode)) then
      md = trim(mode)
    else
      md = 'human'
    end if

    if (md .eq. 'machine') then
      select rank(chains)
        rank(3)
          call store_chains_machine_single(chains, filename)
        rank(4)
          call store_chains_machine_parallel(chains, filename)
      end select
    else if (md .eq. 'human') then
      select rank(chains)
        rank(3)
          call store_chains_human_single(chains, filename)
        rank(4)
          call store_chains_human_parallel(chains, filename)
      end select
    else
      write(*,*) "Optional arugment 'mode' should be either 'human' or 'machine'."
      call exit
    end if

  end subroutine store_chains

  subroutine store_chains_human_single(chains, filename)

    real(wp), intent(in) :: chains(:,:,:)
    character(len=*), intent(in) :: filename

  end subroutine store_chains_human_single

  subroutine store_chains_human_parallel(chains, filename)

    real(wp), intent(in) :: chains(:,:,:,:)
    character(len=*), intent(in) :: filename

  end subroutine store_chains_human_parallel

  subroutine store_chains_machine_single(chains, filename)

    real(wp), intent(in) :: chains(:,:,:)
    character(len=*), intent(in) :: filename

    real(wp), allocatable :: chains2d(:,:)
    integer :: ndim
    integer :: nwalkers
    integer :: nsteps

    integer :: i
    integer :: j
    integer :: k
    integer :: row

    ndim = size(chains, 1)
    nwalkers = size(chains, 2)
    nsteps = size(chains, 3)

    allocate(chains2d(nwalkers * nsteps, ndim))
    chains2d = 0.0e0_wp

    row = 0
    do k = 1, nsteps
       do j = 1, nwalkers
          row = row + 1
          chains2d(row, :) = chains(:, j, k)
       end do
    end do

    open(unit=10, file=filename, form='unformatted', access='stream', status='replace')
      ! Work-around because ifx complains with: write(10) chains2d
      do i = 1, nwalkers * nsteps
        write(10) chains2d(i, :)
      end do
    close(10)

  end subroutine store_chains_machine_single

  subroutine store_chains_machine_parallel(chains, filename)

    real(wp), intent(in) :: chains(:,:,:,:)
    character(len=*), intent(in) :: filename

    real(wp), allocatable :: chains2d(:,:)
    integer :: ndim
    integer :: nwalkers
    integer :: nsteps
    integer :: nthreads

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: row

    nthreads = size(chains, 1)
    ndim = size(chains, 2)
    nwalkers = size(chains, 3)
    nsteps = size(chains, 4)

    allocate(chains2d(nwalkers * nsteps * nthreads, ndim))
    chains2d = 0.0e0_wp

    row = 0
    do k = 1, nsteps
      do j = 1, nwalkers
        do l = 1, nthreads
          row = row + 1
          chains2d(row, :) = chains(l, :, j, k)
        end do
      end do
    end do

    open(unit=10, file=filename, form='unformatted', access='stream', status='replace')
      ! Work-around because ifx complains with: write(10) chains2d
      do i = 1, nwalkers * nsteps * nthreads
        write(10) chains2d(i, :)
      end do
    close(10)

  end subroutine store_chains_machine_parallel

  subroutine store_log_probs(log_probs, filename, mode)

    real(wp), intent(in) :: log_probs(..)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: mode

    character(len=:), allocatable :: md

    if (present(mode)) then
      md = trim(mode)
    else
      md = 'human'
    end if

    if (md .eq. 'machine') then
      select rank(log_probs)
        rank(2)
          call store_log_probs_machine_single(log_probs, filename)
        rank(3)
          call store_log_probs_machine_parallel(log_probs, filename)
      end select
    else if (md .eq. 'human') then
      select rank(log_probs)
        rank(2)
          call store_log_probs_human_single(log_probs, filename)
        rank(3)
          call store_log_probs_human_parallel(log_probs, filename)
      end select
    else
      write(*,*) "Optional arugment 'mode' should be either 'human' or 'machine'."
      call exit
    end if

  end subroutine store_log_probs

  subroutine store_log_probs_human_single(log_probs, filename)

    real(wp), intent(in) :: log_probs(:,:)
    character(len=*), intent(in) :: filename

  end subroutine store_log_probs_human_single

  subroutine store_log_probs_human_parallel(log_probs, filename)

    real(wp), intent(in) :: log_probs(:,:,:)
    character(len=*), intent(in) :: filename

  end subroutine store_log_probs_human_parallel

  subroutine store_log_probs_machine_single(log_probs, filename)

    real(wp), intent(in) :: log_probs(:,:)
    character(len=*), intent(in) :: filename

    real(wp), allocatable :: log_probs_1d(:)
    integer :: ndim
    integer :: nwalkers
    integer :: nsteps

    integer :: i
    integer :: j
    integer :: row

    nwalkers = size(log_probs, 1)
    nsteps = size(log_probs, 2)

    allocate(log_probs_1d(nwalkers * nsteps))
    log_probs_1d = 0.0e0_wp

    row = 0
    do i = 1, nsteps
      do j = 1, nwalkers
        row = row + 1
        log_probs_1d(row) = log_probs(j, i)
      end do
    end do

    open(unit=10, file=filename, form='unformatted', access='stream', status='replace')
      do i = 1, nwalkers * nsteps
        write(10) log_probs_1d(i)
      end do
    close(10)

  end subroutine store_log_probs_machine_single

  subroutine store_log_probs_machine_parallel(log_probs, filename)

    real(wp), intent(in) :: log_probs(:,:,:)
    character(len=*), intent(in) :: filename

    real(wp), allocatable :: log_probs_1d(:)
    integer :: ndim
    integer :: nwalkers
    integer :: nsteps
    integer :: nthreads

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: row

    nthreads = size(log_probs, 1)
    nwalkers = size(log_probs, 2)
    nsteps = size(log_probs, 3)

    allocate(log_probs_1d(nwalkers * nsteps * nthreads))
    log_probs_1d = 0.0e0_wp

    row = 0
    do k = 1, nsteps
      do j = 1, nwalkers
        do l = 1, nthreads
          row = row + 1
          log_probs_1d(row) = log_probs(l, j, k)
        end do
      end do
    end do

    open(unit=10, file=filename, form='unformatted', access='stream', status='replace')
      do i = 1, nwalkers * nsteps * nthreads
        write(10) log_probs_1d(i)
      end do
    close(10)

  end subroutine store_log_probs_machine_parallel

end module copa__store
