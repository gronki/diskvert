module grid

  use iso_fortran_env, only: r64 => real64
  implicit none

contains

  ! linear grid
  elemental function space_linear(i,n,h) result(f)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: h
    real(r64) :: f,x
    x = real(i - 1, r64) / real(n - 1, r64)
    f = x * h
  end function

  ! logarithm-like grid (starting at 0)
  ! goes from 0 to h and it's more linear for h < 1
  ! and more logarithmic for h > 1.
  elemental function space_linlog(i,n,h) result(f)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: h
    real(r64) :: f,x
    x = real(i - 1, r64) / real(n - 1, r64)
    f = (1 + h) ** x - 1
  end function

  ! logarithm-like grid (starting at 0)
  ! goes from 0 to h and it's more linear for h < 1
  ! and more logarithmic for h > 1.
  elemental function space_asinh(i,n,h) result(f)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: h
    real(r64) :: f,x
    x = real(i - 1, r64) / real(n - 1, r64)
    f = sinh( x * asinh(h) )
  end function

  elemental function space_pow2(i,n,k) result(f)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: k
    real(r64) :: f,t
    t = real(i - 1, r64) / real(n - 1, r64)
    f = (2 * t + (k - 1) * t**2) / (k + 1)
  end function

  ! linear-logarithmic log
  ! linear steps from 0 to y1, and logarithmic from y1 to y2
  ! it is assumed that y1 = 1 and y2 > 1
  pure subroutine space_linlog2(x, y2)
    real(r64), intent(out) :: x(:)
    real(r64), intent(in) :: y2
    real(r64) :: y1, n2f
    integer :: n1, n2, n, i
    
    n = size(x)
    n2f = n - 1
    y1 = 1.0

    do i = 1, 50
      N2f = N - 1 / (1 - (Y1 / Y2)**(1 / N2f))
    end do

    n2 = min(floor(n2f), n - 1)
    n1 = n - n2

    do i = 1, n1
      x(i) = (i - 1) * (y1 / n1)
    end do

    do i = n1, n
      x(i) = y1 * (y2 / y1)**((i - n1) / real(n2, r64))
    end do
  end subroutine

end module