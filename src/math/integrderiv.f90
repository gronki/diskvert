module slf_integrate

  use iso_fortran_env, only: r32 => real32, r64 => real64
  implicit none

  interface integrate
    module procedure integrate32
    module procedure integrate64
  end interface integrate

contains

  pure real(r32) function integrate32(y,x) result(intg)
    real(r32), dimension(:), intent(in) :: x, y
    integer :: n

    if (size(x) /= size(y)) error stop "size(x) /= size(y)"

    n = size(x)
    intg = sum((y(2:n) + y(1:n-1)) * (x(2:n) - x(1:n-1))) / 2
  end function

  pure real(r64) function integrate64(y,x) result(intg)
    real(r64), dimension(:), intent(in) :: x, y
    integer :: n

    if (size(x) /= size(y)) error stop "size(x) /= size(y)"

    n = size(x)
    intg = sum((y(2:n) + y(1:n-1)) * (x(2:n) - x(1:n-1))) / 2
  end function

end module slf_integrate

!------------------------------------------------------------------------------!

module slf_deriv

  use iso_fortran_env, only: r64 => real64
  implicit none

contains

  subroutine deriv(x,yin,yout)
    real(r64), intent(in) :: x(:), yin(:)
    real(r64), intent(out) :: yout(:)
    real(r64) :: dx2
    integer :: i,n

    n = size(yout)

    if ( n .lt. 3 ) then
      write (0,*) 'n should be at least 3!'
      return
    end if

    do i=2,n-1
      yout(i) = (yin(i+1)-yin(i-1))/(x(i+1)-x(i-1))
    end do
    dx2 = 2 * ( yin(3) + yin(1) - 2 * yin(2) ) / ( (x(3)-x(2))**2 + (x(2)-x(1))**2 )
    yout(1) = yout(2) + (x(1)-x(2))*dx2
    dx2 = 2 * ( yin(n) + yin(n-2) - 2 * yin(n-1) ) / ( (x(n)-x(n-1))**2 + (x(n-1)-x(n-2))**2 )
    yout(n) = yout(n-1) + (x(n)-x(n-1))*dx2
  end subroutine

  subroutine deriv2(x,yin,yout)
    real(r64), intent(in) :: x(:), yin(:)
    real(r64), intent(out) :: yout(:)
    integer :: i,n
    n = size(yout)

    do i=2,n-1
      yout(i) = 2 * ( yin(i+1) + yin(i-1) - 2 * yin(i) ) / ( (x(i+1)-x(i))**2 + (x(i-1)-x(i))**2 )
    end do
    yout(1) = yout(2) + ( x(1) - x(2) )*( yout(3)-yout(2) )/( x(3)-x(2) )
    yout(n) = yout(n-1) + ( x(n) - x(n-1) )*( yout(n-1)-yout(n-2) )/( x(n-1)-x(n-2) )
  end subroutine

   !----------------------------------------------------------------------------!
  ! computes log gradient of any function

  pure subroutine loggrad(x, y, d)
    real(r64), intent(in) :: x(:), y(:)
    real(r64), intent(out) :: d(:)
    integer :: i

    if (size(x) /= size(y) .or. size(d) /= size(y)) error stop

    do concurrent (i = 2:size(y)-1)
      d(i) = (y(i+1) - y(i-1)) * (x(i+1) + x(i-1)) &
      &   / ((y(i+1) + y(i-1)) * (x(i+1) - x(i-1)) )
    end do

    d(1) = d(2)
    d(size(d)) = d(size(d) - 1)
  end subroutine


end module
