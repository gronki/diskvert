module alphadisk

    use iso_fortran_env, only: r64 => real64
    use iso_c_binding

    use globals
    use slf_cgs
    use grid
    use slf_rk4integr
    use slf_eulerintegr
    use slf_threshold, only: thrtanh
    use rk4settings
    use fileunits
    use heatbalance
    use ss73solution, only: ss73_estim, ss73_refin

    implicit none

    real(r64), private :: radius, alpha
    real(r64), private :: omega, Facc, Teff, zscale
    real(r64), private :: rhoc, tempc, hdisk

    integer, parameter :: ny = 4,   &
        &   c_Pgas = 1, c_Prad = 2, c_Frad = 3, c_tau = 4
    integer, parameter :: na = 10,   &
        &   c_zscal   = 1, &
        &   c_Tgas    = 2, &
        &   c_rho     = 3, &
        &   c_Trad    = 4, &
        &   c_heat    = 5, &
        &   c_heat_max = 6, &
        &   c_compW   = 7, &
        &   c_compY   = 8, &
        &   c_kabs = 9, c_ksct = 10

contains

    subroutine init_alpha(mb,md,r,alph) bind(C)
        real(r64), intent(in) :: mb,md,r,alph
        mbh = mb
        mdot = md
        radius = r
        alpha = alph
        call cylinder(mbh, mdot, radius, rschw, omega, Facc, Teff, zscale)
        call ss73_estim(mbh, mdot, radius, alpha, rhoc, tempc, hdisk)
        call ss73_refin(mbh, mdot, radius, alpha, rhoc, tempc, hdisk)
    end subroutine

    SUBROUTINE ALF_GETN(NY1,NA1) BIND(C)
      INTEGER(C_INT), INTENT(OUT) :: NY1,NA1
      NY1 = NY
      NA1 = NA
    END SUBROUTINE

    subroutine run_alpha(z,nz,y,dy,a,nmax) bind(C)
        ! ilosc przedzialow obliczeniowych
        integer, intent(in) :: nz
        real(r64), intent(inout), dimension(nz) :: z
        real(r64), intent(inout), dimension(ny,nz) :: y,dy
        real(r64), intent(inout), dimension(na,nz) :: a
        integer(c_int), intent(out) :: nmax
        real(r64) :: h_hi, h_lo, h
        integer :: i,it,niter

        h = 4 * hdisk / zscale
        niter = ceiling((log10(1e1-1e-1)-log10(max_iteration_error))/log10(2.0))
        h_hi = h * 10
        h_lo = h / 10

        write (uerr,'(A10,A10,A12,A7,A5)') "ITER", "H", "FLX", "MAX", "ACTN"

        do it=1,niter
            h = sqrt(h_hi*h_lo)

            forall (i = 1:nz)  z(i) = space_linlog_r(i,nz,h) * zscale

            y(c_Frad,1) = Facc
            y(c_Prad,1) = Facc * 2 / (3 * cgs_c)
            y(c_Pgas,1) = Teff * 1e-16 * cgs_k_over_mh
            y(c_tau,1)  = epsilon(1d0)

            if ( cfg_euler_integration ) then
                call eulerintegr(z, y(:,1), f, na, y, dy, a, nmax)
            else
                call rk4integr(z, y(:,1), f, na, y, dy, a, nmax)
            end if

            if ( nmax < nz ) then
                h_hi = h
                write (uerr,'(I10,F10.3,Es12.4,F7.1,A5)') it, h, y(c_Frad,nmax) / Facc, 100.0 * nmax / nz, '\/'
            else if (y(c_Frad,nz) > 1e-10 * Facc) then
                h_lo = h
                write (uerr,'(I10,F10.3,Es12.4,F7.1,A5)') it, h, y(c_Frad,nmax) / Facc, 100.0 * nmax / nz, '/\'
            else
                write (uerr,'("iteration finished")')
                exit
            end if
        end do
    end subroutine


    subroutine f(z,y,dy,a,abort)

        real(r64), intent(in) :: z, y(:)
        real(r64), intent(inout) :: dy(size(y)), a(:)
        logical, intent(inout) :: abort

        real(r64) :: compsw
        real(r64), parameter :: toler = 1
        real(r64) ::  epsi, taues

        if ( y(c_Frad) < 0 ) then
            abort = .TRUE.
        end if

        a(c_zscal)  = z / zscale
        a(c_Trad)   = ( (3 * cgs_c) / (4 * cgs_stef) * y(c_Prad)) ** 0.25d0
        a(c_heat)   = alpha * omega * y(c_Pgas)

        if ( cfg_temperature_method == EQUATION_DIFFUSION ) then
          a(c_heat) = a(c_heat) + alpha * omega * y(c_prad)
        end if

        a(c_heat_max) = 12 * cgs_kapes * y(c_Prad)  &
                & * ( y(c_Pgas) * miu * cgs_mhydr ) &
                & / ( cgs_mel * cgs_c )

        a(c_compW) = a(c_heat) / ( a(c_heat_max) - a(c_heat) )
        if ( a(c_compW) < 0 )   a(c_compW) = 0

        select case (cfg_temperature_method)
        case (EQUATION_DIFFUSION)
            a(c_Tgas)   = a(c_Trad)
        case (EQUATION_COMPTON)
            a(c_Tgas) = a(c_Trad) * (1 + a(c_compW))
            a(c_rho) = y(c_Pgas) * miu / ( cgs_k_over_mh * a(c_Tgas) )

            epsi = 1 / (1 + fkesp(a(c_rho), a(c_Tgas)) / fkabp(a(c_rho), a(c_Tgas)))
            taues = (1. - epsi)/sqrt(epsi)
            a(c_compY) = 4 * cgs_boltz * a(c_tgas) / ( cgs_mel * cgs_c**2 )  &
                & * max( taues, taues**2 )
            compsw = thrtanh( (log(a(c_compY)) - log(1.0))/log(1+toler) )

            a(c_Tgas) = a(c_trad) * ( 1 + compsw * a(c_compW) )
        case (EQUATION_BALANCE)
            a(c_Tgas) = a(c_trad)
            call heatbil(a(c_tgas), a(c_trad), y(c_pgas), a(c_heat))
        case default
            error stop "incorrect value of cfg_temperature_method"
        end select

        a(c_rho) = y(c_Pgas) * miu / ( cgs_k_over_mh * a(c_Tgas) )

        epsi = 1 / (1 + fkesp(a(c_rho),a(c_Tgas)) / fkabs(a(c_rho),a(c_Tgas)))
        taues = (1 - epsi) / sqrt(epsi)
        a(c_compY) = 4 * cgs_boltz * a(c_tgas) &
            & / ( cgs_mel * cgs_c**2 )  &
            & * max( taues, taues**2 )

        a(c_kabs) = fkabs(a(c_rho),a(c_Tgas))
        a(c_ksct) = fkesp(a(c_rho),a(c_Tgas))

        dy(c_Frad)  = a(c_heat)
        dy(c_Prad)  = - cgs_kapes * a(c_rho) / cgs_c * y(c_Frad)
        dy(c_Pgas)  = - omega**2 * a(c_rho) * z - dy(c_Prad)
        dy(c_tau)   = - a(c_kabs) * a(c_rho)

    end subroutine

end module
