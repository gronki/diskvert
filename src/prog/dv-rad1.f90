program dv_rad1

    use iso_fortran_env, only: r64 => real64
    use slf_cgs
    use ss73solution
    use globals
    use alphadisk
    use fileunits
    use confort
    use settings

    implicit none

    integer :: i, j
    integer, parameter :: n = 2**9
    real(r64) :: r12, r23, alpha
    real(r64), dimension(n) :: r
    real(r64), dimension(n) :: rho1, T1, H1
    real(r64), dimension(n) :: rho2, T2, H2
    real(r64), dimension(n) :: rho3, T3, H3
    type(config) :: cfg

    integer, parameter :: nz = 2**11
    real(r64) :: z(nz), y(ny,nz), dy(ny,nz), ap(na,nz)
    integer :: nmax


    call rdargvgl()
    call mincf_read(cfg)
    call rdconfgl(cfg)
    call rdconf(cfg)
    call mincf_free(cfg)

    call ss73_zonebounds(mbh,mdot,alpha,r12,r23)

    do i = 1, n
        r(i) = 3 * exp( i / real(n) * (log(r23*10)-log(3d0)) )
    end do

    call ss73_sel(mbh, mdot, r, alpha, r12, r23, rho1, T1, H1)

    rho2 = rho1
    T2 = T1
    H2 = H1

    compute_precise_model: do i = 1,n
        call ss73_refin(mbh, mdot, r(i), alpha, rho2(i), T2(i), H2(i))

        call init_alpha(mbh, mdot, r(i), alpha)
        call run_alpha(z,nz,y,dy,ap,nmax)
        rho3(i) = ap(c_rho,nz)
        T3(i) = ap(c_Tgas,nz)
        H3(i) = z(1) / 2
        search_height: do j = 1,nz
            if (ap(c_rho,j) .ge. rho3(i) / 2) then
                H3(i) = z(j)
                exit search_height
            end if
        end do search_height
    end do compute_precise_model

    do i = 1,n
        write (6, '(I7,10(1X,Es11.4))') i, r(i), rho1(i), T1(i), H1(i), &
            rho2(i), T2(i), H2(i), rho3(i), T3(i), H3(i)
    end do

contains

    subroutine rdconf(cfg)
        type(config), intent(inout) :: cfg
        character(128) :: buf

        if ( .not. mincf_exists(cfg,'alpha') ) then
            error stop "alpha key is required!"
        end if

        call mincf_get(cfg,'alpha',buf)
        read (buf,*) alpha
    end subroutine

end program
