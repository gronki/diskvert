!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module slf_findzer

    use iso_fortran_env, only: r64 => real64
    implicit none

    interface findzer
        module procedure findzer_single, findzer_array
    end interface

    integer, parameter :: findzer_iter = 64
    private :: ramp

contains

    subroutine findzer_single(x, xlo, xhi, delx, f)

        real(r64), intent(inout) :: x
        real(r64), intent(in) :: xlo,xhi,delx

        interface
            pure subroutine f(x,y,dy)
                import r64
                real(r64), intent(in) :: x
                real(r64), intent(out) :: y
                real(r64), intent(out), optional :: dy
            end subroutine
        end interface

        integer :: i
        real(r64) :: dx,yx,y,x1

        main_loop : do i = 1,FINDZER_ITER

            call f(x,y,yx)
            if (yx /= 0) then
                dx = - y / yx
            elseif ( dx == 0 ) then
                exit main_loop
            end if

            x1 = x + dx * ramp(i,findzer_iter)
            if ( x1 > xhi ) x1 = (x + xhi) / 2
            if ( x1 < xlo ) x1 = (x + xlo) / 2
            x = x1

            if (abs(dx) < abs(delx)) exit main_loop

        end do main_loop

    end subroutine

    subroutine findzer_array(x, nx, xlo, xhi, delx, f)

        integer, intent(in) :: nx
        real(r64), intent(inout), dimension(nx) :: x
        real(r64), intent(in), dimension(nx)  :: xlo,xhi
        real(r64), intent(in)  :: delx

        interface
            pure subroutine f(x,y,dy)
                import r64
                real(r64), intent(in), dimension(:)  :: x
                real(r64), intent(out), dimension(size(x))  :: y
                real(r64), intent(out), dimension(size(x)) , optional :: dy
            end subroutine
        end interface

        integer :: i
        real(r64), dimension(nx) :: dx,yx,y,x1
        logical, dimension(nx) :: converged

        converged = .false.
        dx = 0

        main_loop : do i = 1,FINDZER_ITER

            call f(x,y,yx)

            where ( .not. converged )

                where ( yx /= 0 )
                    dx = - y / yx
                elsewhere ( dx == 0 )
                    converged = .true.
                end where

                x1 = x + dx * ramp(i,findzer_iter)
                where ( x1 > xhi ) x1 = (x + xhi) / 2
                where ( x1 < xlo ) x1 = (x + xlo) / 2
                x = x1

                converged = abs(dx) < abs(delx)
            end where

            if (all(converged)) exit main_loop

        end do main_loop

    end subroutine

    elemental function ramp(i,n) result(y)
      integer, intent(in) :: i,n
      real(r64) :: t,y
      t = merge(real(i) / n, 1.0, i .lt. n)
      y = (3 - 2*t) * t**2
    end function

end module

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module slf_findzermulti

    use ieee_arithmetic
    use iso_fortran_env, only: r64 => real64
    implicit none

contains

    subroutine findzer_multi(x, nx, xlo, xhi, delx, nb, f)

        real(r64), intent(inout) :: x(:)
        real(r64), intent(in) :: xlo,xhi,delx
        integer, intent(in) :: nb
        integer, intent(out) :: nx

        interface
            pure subroutine f(x,y,dy)
                import r64
                real(r64), intent(in) :: x
                real(r64), intent(out) :: y
                real(r64), intent(out), optional :: dy
            end subroutine
        end interface

        real(r64) :: y(nb), dy(nb), x0(nb)
        integer :: i


        do concurrent (i=1:nb)
            x0(i) = (i-1)/(nb-1.) * (xhi - xlo) + xlo
        end do

        call f(x0(1),y(1))
        nx = 0
        scan_for_solutions: do i = 2, nb
            call f(x0(i),y(i))
            if ( y(i-1) * y(i) < 0 ) then
                call linsect(x(nx+1), x0(i-1), x0(i), delx)
                nx = nx + 1
            end if
            if ( nx >= size(x) ) exit scan_for_solutions
        end do scan_for_solutions

    contains

        subroutine linsect(x,xlo0,xhi0,delx)
            real(r64), intent(in) :: xhi0,xlo0,delx
            real(r64), intent(out) :: x
            real(r64) :: xhi,xlo,yhi,ylo,y,s
            integer :: i,n

            n = log(abs((xhi0-xlo0)/delx))/log(2.0)

            xlo = xlo0
            xhi = xhi0

            sect_loop: do i = 1, n
                call f(xlo,ylo)
                call f(xhi,yhi)
                s = sign(real(1,kind=r64),yhi-ylo)

                x = ( yhi * xlo - ylo * xhi ) / ( yhi - ylo )
                call f(x,y)

                if ( y * s > 0 ) then
                    xhi = x
                elseif ( y * s < 0 ) then
                    xlo = x
                else
                    exit sect_loop
                end if

                if ( abs(xhi-xlo) < abs(delx) ) exit sect_loop
            end do sect_loop

        end subroutine

    end subroutine

end module

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module slf_bisect

    use iso_fortran_env, only: r64 => real64
    implicit none

contains

    subroutine bisect(f,x,xlo0,xhi0,delx)
        real(r64), intent(in) :: xhi0,xlo0,delx
        real(r64), intent(out) :: x
        real(r64) :: xhi,xlo,yhi,ylo,y,s
        integer :: i,n

        interface
            pure subroutine f(x,y)
                import r64
                real(r64), intent(in) :: x
                real(r64), intent(out) :: y
            end subroutine
        end interface

        n = log(abs((xhi0-xlo0)/delx))/log(2.0)

        xlo = xlo0
        xhi = xhi0
        call f(xlo,ylo)
        call f(xhi,yhi)
        s = sign(real(1,r64),yhi-ylo)

        sect_loop: do i = 1, n

            x = ( xhi + xlo ) / 2
            call f(x,y)

            if ( y * s > 0 ) then
                xhi = x
            elseif ( y * s < 0 ) then
                xlo = x
            else
                exit sect_loop
            end if

            if ( abs(xhi-xlo) < abs(delx) ) exit sect_loop
        end do sect_loop

    end subroutine

end module

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module slf_linsect

    use iso_fortran_env, only: r64 => real64
    implicit none

contains

    subroutine linsect(f,x,xlo0,xhi0,delx)
        real(r64), intent(in) :: xhi0,xlo0,delx
        real(r64), intent(out) :: x
        real(r64) :: xhi,xlo,yhi,ylo,y,s
        integer :: i,n

        interface
            pure subroutine f(x,y)
                import r64
                real(r64), intent(in) :: x
                real(r64), intent(out) :: y
            end subroutine
        end interface

        n = log(abs((xhi0-xlo0)/delx))/log(2.0)

        xlo = xlo0
        xhi = xhi0

        sect_loop: do i = 1, n
            call f(xlo,ylo)
            call f(xhi,yhi)
            s = sign(real(1,r64),yhi-ylo)

            x = ( yhi * xlo - ylo * xhi ) / ( yhi - ylo )
            call f(x,y)

            if ( y * s > 0 ) then
                xhi = x
            elseif ( y * s < 0 ) then
                xlo = x
            else
                exit sect_loop
            end if

            if ( abs(xhi-xlo) < abs(delx) ) exit sect_loop
        end do sect_loop

    end subroutine

end module

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
