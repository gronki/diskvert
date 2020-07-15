module ranges

  use globals
  implicit none

contains

  pure subroutine linrange(x, xlo, xhi)
    real(r64), intent(out) :: x(:)
    real(r64), intent(in) :: xlo, xhi
    integer :: n, i
    real(r64) :: t

    n = size(x)

    do concurrent (i = 1:n)
      t = real(i - 1, r64) / real(n - 1, r64)
      x(i) = xlo * (1 - t) + xhi * t
    end do
  end subroutine

  pure subroutine logrange(x, xlo, xhi)
    real(r64), intent(out) :: x(:)
    real(r64), intent(in) :: xlo, xhi

    call linrange(x, log(xlo), log(xhi))
    x(:) = exp(x)
  end subroutine

end module