module slf_interpol

  use iso_fortran_env, only: r64 => real64
  implicit none

contains

  pure subroutine interpol(x_in, y_in, x_out, y_out)
    real(r64), dimension(:), intent(in) :: x_in, y_in
    real(r64), intent(in) :: x_out
    real(r64), intent(out) :: y_out
    real(r64) :: x0,x1,y0,y1,t,control
    integer :: i, i_nearest, n

    if (size(y_in) /= size(x_in)) error stop "size(y_in) /= size(x_in)"

    n = size(y_in)
    t = abs(x_in(1)-x_in(n))

    ! find the nearest point
    control = x_in(2) - x_in(1)
    do i=1,n
      if (i .gt. 1) then
        if ( control * (x_in(i)-x_in(i-1)) .le. 0 ) then
          error stop 'interpol: x-axis points should be either in ascending or descending order'
        endif
      endif
      if ( (abs(x_in(i)-x_out) .lt. abs(t)) .or. (i .eq. 1) ) then
        i_nearest = i
        t = x_in(i)-x_out
      end if
    end do

    if ( t .eq. 0. ) then
      y_out = y_in(i_nearest)
      return
    else if (t .gt. 0) then
      if ( i_nearest .eq. 1 ) then
        y_out = y_in(1)
        return
      endif
      x1 = x_in(i_nearest)
      y1 = y_in(i_nearest)
      x0 = x_in(i_nearest-1)
      y0 = y_in(i_nearest-1)
    else
      if ( i_nearest .eq. n ) then
        y_out = y_in(n)
        return
      endif
      x0 = x_in(i_nearest)
      y0 = y_in(i_nearest)
      x1 = x_in(i_nearest+1)
      y1 = y_in(i_nearest+1)
    end if

    t = (x1 - x_out) / (x1-x0)
    y_out = t * y0 + (1.-t) * y1

  end subroutine

  !----------------------------------------------------------------------------!

  pure function interpolf(x,y,x0) result(y0)
    real(r64), intent(in) :: x(:), y(:), x0
    real(r64) :: y0
    call interpol(x,y,x0,y0)
  end function

end module
