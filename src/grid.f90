module grid

  use iso_fortran_env, only: r64 => real64
  implicit none

contains

  elemental function space_linear(i,n,h) result(f)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: h
    real(r64) :: f,x
    x = real(i - 1, r64) / real(n - 1, r64)
    f = x * h
  end function

  elemental function space_linlog(i,n,h) result(f)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: h
    real(r64) :: f,x
    x = real(i - 1, r64) / real(n - 1, r64)
    f = (1 + h) ** x - 1
  end function

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

end module
