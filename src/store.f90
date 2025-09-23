module copa__store

  use copa__config, only : wp

  implicit none

  private :: store_chains_human
  private :: store_chains_machine

  public :: store_chains

contains

  subroutine store_chains(chains, filename, mode)

    real(wp), intent(in) :: chains(:,:,:)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: mode

    character(len=:), allocatable :: md

    if (present(mode)) then
      md = trim(mode)
    else
      md = 'human'
    end if

    select case(md)
      case ('human')
        call store_chains_human(chains, filename)
      case ('machine')
        call store_chains_machine(chains, filename)
      case default
        write(*,*) "Optional arugment 'mode' should be either 'human' or 'machine'."
        call exit
    end select

  end subroutine store_chains

  subroutine store_chains_human(chains, filename)

    real(wp), intent(in) :: chains(:,:,:)
    character(len=*), intent(in) :: filename

  end subroutine store_chains_human

  subroutine store_chains_machine(chains, filename)

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
          chains2d(row,:) = chains(:,j, k)
       end do
    end do

    open(unit=10, file=filename, form='unformatted', access='stream', status='replace')
      ! Work-around because ifx complains with: write(10) chains2d
      do i = 1, nwalkers * nsteps
        write(10) chains2d(i,:)
      end do
    close(10)

  end subroutine store_chains_machine

end module copa__store
