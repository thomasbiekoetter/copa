module copa__prior_functions

  use copa__config, only : wp

  implicit none

  private

  real(wp), parameter :: huge_negative = -1.0e45_wp
  real(wp), parameter :: gaussian_f = 1.0e0_wp
  real(wp), parameter :: gaussian_eps = 1.0e-10_wp

  public :: uniform_prior
  public :: log_uniform_prior
  public :: gaussian_prior

contains

  subroutine uniform_prior(theta, lower, upper, logp)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(in) :: lower(:)
    real(wp), intent(in) :: upper(:)
    real(wp), intent(out) :: logp

    if (all(theta >= lower) .and. all(theta < upper)) then
      logp = 0.0e0_wp
    else
      logp = huge_negative
    end if

  end subroutine uniform_prior

  subroutine log_uniform_prior(theta, lower, upper, logp)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(in) :: lower(:)
    real(wp), intent(in) :: upper(:)
    real(wp), intent(out) :: logp

    integer :: i

    ! Log-uniform prior (Jeffrey-like)
    if (any(lower <= 0.0e0_wp)) then
      write(*,*) "Lower limit of thetas must be positive for log-uniform prior."
      call exit
    end if
    logp = huge_negative
    if (all((theta >= lower) .and. (theta <= upper))) then
      logp = 0.0e0_wp
      do i = 1, size(theta)
        ! Metropolis–Hastings acceptance rule -> ln (not log10)
        logp = logp - log(theta(i)) - (log(log(upper(i)) -  &
          log(lower(i))))
      end do
    end if

  end subroutine log_uniform_prior

  subroutine gaussian_prior(theta, mu, logp, sigma, f_sigma)

    real(wp), intent(in) :: theta(:)
    real(wp), intent(in) :: mu(:)
    real(wp), intent(out) :: logp
    real(wp), intent(in), optional :: sigma(:)
    real(wp), intent(in), optional :: f_sigma(:)

    real(wp), allocatable :: sig(:)
    real(wp), allocatable :: f(:)

    allocate(sig(size(theta)))
    allocate(f(size(theta)))

    if (present(sigma)) then
      sig = sigma
    else
      if (present(f_sigma)) then
        f = f_sigma
      else
        f = gaussian_f
      end if
      sig = max(f * abs(mu), gaussian_eps)
    end if

    logp = -0.5_wp * sum(((theta - mu) / sig) ** 2)

  end subroutine gaussian_prior

end module copa__prior_functions
