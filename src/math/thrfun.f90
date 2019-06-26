module slf_threshold

  use iso_fortran_env, only: r64 => real64
  implicit none

  real(r64), parameter, private :: pi = 4*atan(real(1,r64))

contains

  elemental real(r64) function THRSTEP(x) result(y)
    real(r64), intent(in) :: x
    if ( x .ge. 0 ) then
      y = 1.
    else
      y = 0.
    end if
  end function

  elemental real(r64) function THRLIN(x) result(y)
    real(r64), intent(in) :: x
    if ( x < -0.5_r64) then
      y = 0
    elseif ( x > 0.5_r64) then
      y = 1
    else
      y = x + 0.5_r64
    end if
  end function

  elemental real(r64) function THR2POLY(x,ord) result(y)
    real(r64), intent(in) :: x
    integer, intent(in), optional :: ord
    integer :: n
    real(r64) :: a, b

    n = 3
    if ( present(ord) )   n = ord

    a = THRLIN(0.5 + 2*x/n)
    b = THRLIN(0.5 - 2*x/n)

    y = (a**n - b**n + 1) / 2
  end function

  elemental real(r64) function THR4POLY(x,ord) result(y)
    real(r64), intent(in) :: x
    integer, intent(in), optional :: ord
    integer :: n
    real(r64) :: a, b, c

    n = 3
    if ( present(ord) )   n = ord

    a = THRLIN((3 + x*4)/2)
    b = 2*THRLIN(x) - 1
    c = THRLIN((3 - x*4)/2)

    y = (a**n - c**n + b*(2*n - abs(b)**(n-1)) + 2*n)/(4*n)

  end function

  elemental real(r64) function THRSIN(x) result(y)
    real(r64), intent(in) :: x
    real(r64) :: z

    z = 2*THRLIN(x / 2) - 1
    y =  (sin( pi * z) / pi + z + 1) / 2
  end function

  elemental real(r64) function THRATAN(x) result(y)
    real(r64), intent(in) :: x
    y = atan(x*pi) / pi + 0.5
  end function

  elemental real(r64) function thrsqrt(x) result(y)
    real(r64), intent(in) :: x
    y = 0.5_r64 +  x / sqrt( 1 + 4 * x**2 )
  end function

  elemental real(r64) function thrtanh(x) result(y)
    real(r64), intent(in) :: x
    y = 1d0 / ( 1d0 + exp(-4*x) )
  end function

end module

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module slf_ramps

  use iso_fortran_env, only: r64 => real64
  implicit none

contains

  !----------------------------------------------------------------------------!
  ! this is a linear ramp

  elemental function ramp1(i,n) result(y)
    integer, intent(in) :: i,n
    real(r64) :: y
    y = real(max(min(i, n), 0), r64) / n
  end function

  !----------------------------------------------------------------------------!
  ! this is a very steep ramp, for later iterations

  elemental function ramp2(i,n) result(y)
    integer, intent(in) :: i,n
    real(r64) :: t,y
    t = real(max(min(i, n), 0), r64) / n
    y = t * (2 - t)
  end function

  !----------------------------------------------------------------------------!
  ! this is a very smooth, s-shaped ramp
  ! slow convergence but good if we are away from the solution

  elemental function ramp3(i,n) result(y)
    integer, intent(in) :: i,n
    real(r64) :: t,y
    t = real(max(min(i, n), 0), r64) / n
    y = (3 - 2*t) * t**2
  end function

  !----------------------------------------------------------------------------!
  ! this is a very smooth ramp where you can control the starting
  ! level as well as slope of the starting point. It may be applied
  ! for very unstable equations (A = 0, B = 0) or for a very quick
  ! convergence (A = 0.25, B = 1)

  elemental real(r64) function ramp4r(x, A, B) result(y)
    real(r64), intent(in) :: x
    real(r64), intent(in) :: A,B
    y = A + (1 - A) * B * x           &
    + 3 * (A - 1) * (B - 2) * x**2  &
    - (A - 1) * (3*B - 8) * x**3    &
    + (A - 1) * (B - 3) * x**4
  end function

  elemental function ramp4(i,n,A,B) result(y)
    integer, intent(in) :: i,n
    real(r64), intent(in) :: A,B
    real(r64) :: x,y
    x = real(max(min(i, n), 0), r64) / n
    y = ramp4r(x, A, B)
  end function

  !----------------------------------------------------------------------------!
  ! this ramp starts smoothly, rises rapidly, and then slowly saturates to 1

  elemental real(r64) function ramp5r(x) result(y)
    real(r64), intent(in) :: x
    y = (6 + x * (x * 3 - 8)) * x**2
  end function

  elemental function ramp5(i,n) result(y)
    integer, intent(in) :: i,n
    real(r64) :: t,y
    t = real(max(min(i, n), 0), r64) / n
    y = ramp5r(t)
  end function

  !----------------------------------------------------------------------------!

  elemental real(r64) function ramp6r(x) result(y)
    real(r64), intent(in) :: x
    y = x**2 * (5 + x * (10 + x * (x * 8 - 20))) / 3
  end function

  elemental function ramp6(i,n) result(y)
    integer, intent(in) :: i,n
    real(r64) :: t,y
    t = real(max(min(i, n), 0), r64) / n
    y = ramp6r(t)
  end function

  !----------------------------------------------------------------------------!

end module slf_ramps
