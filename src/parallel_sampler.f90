module copa__parallel_sampler

  use copa__config, only : wp
  use copa__config, only : parallel_method_default
  use copa__independent_parallel_sampler, only : run_independent_parallel_sampler
  use copa__redblack_parallel_sampler, only : run_redblack_parallel_sampler

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
    ndim, log_prior, log_like,  &
    method,  &
    nwalkers, nsteps, nthreads, nensembles,  &
    ranges, walkers, chains, log_probs)
    integer, intent(in) :: ndim
    procedure(log_prior_abstract) :: log_prior
    procedure(log_like_abstract) :: log_like
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: nwalkers
    integer, intent(in), optional :: nsteps
    integer, intent(in), optional :: nthreads
    integer, intent(in), optional :: nensembles
    real(wp), intent(in), optional :: ranges(:,:)
    real(wp), intent(out), allocatable, optional :: walkers(:,:,:)
    real(wp), intent(out), allocatable, optional :: chains(:,:,:,:)
    real(wp), intent(out), allocatable, optional :: log_probs(:,:,:)

    character(len=:), allocatable :: method_

    method_ = parallel_method_default
    if (present(method)) method_ = method

    select case (trim(adjustl(method_)))
      case ('independent')
        call run_independent_parallel_sampler( &
          ndim, log_prior, log_like, nwalkers, nsteps, nthreads, &
          ranges, walkers, chains, log_probs)
      case ('redblack')
        call run_redblack_parallel_sampler( &
          ndim, log_prior, log_like, nwalkers, nsteps, nthreads, nensembles, &
          ranges, walkers, chains, log_probs)
      case default
        write(*,'(3a)') "Unknown method '", trim(method_), &
          "'. Valid options: 'independent', 'redblack'."
        call exit
    end select

  end subroutine run_parallel_sampler

end module copa__parallel_sampler
