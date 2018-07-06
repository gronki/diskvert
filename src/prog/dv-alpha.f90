program dv_alpha

    use iso_fortran_env, only: r64 => real64
    use confort
    use globals
    use results
    use settings
    use alphadisk
    use fileunits

    implicit none

    character(len=12), dimension(ny) :: y_labels
    character(len=12), dimension(na) :: a_labels
    integer :: nmax
    real(r64), allocatable, dimension(:) :: z
    real(r64), allocatable, dimension(:,:) :: y,dy,a
    real(r64) :: radius, alpha
    type(config) :: cfg
    integer, parameter :: upar = 92


    y_labels(c_Pgas)    = 'Pgas'
    y_labels(c_Prad)    = 'Prad'
    y_labels(c_Frad)    = 'Frad'
    y_labels(c_tau)     = 'tau_es'
    a_labels(c_rho)     = 'rho'
    a_labels(c_Tgas)    = 'Tgas'
    a_labels(c_Trad)    = 'Trad'
    a_labels(c_zscal)   = 'h'
    a_labels(c_heat)    = 'Gamma'
    a_labels(c_heat_max)= 'Gammamax'
    a_labels(c_compW)   = 'compW'
    a_labels(c_compY)   = 'compY'
    a_labels(c_kabs)    = 'kabs'
    a_labels(c_ksct)    = 'ksct'

    ngrid = 4096

    call rdargvgl
    call rdargvrk4
    call mincf_read(cfg)
    call rdconfgl(cfg)
    call rdconf(cfg)
    call mincf_free(cfg)

    allocate(z(ngrid))
    allocate(y(ny,ngrid), dy(ny,ngrid), a(na,ngrid))

    call init_alpha(mbh, mdot, radius, alpha)
    call run_alpha(z,ngrid,y,dy,a,nmax)

    open(upar, file = trim(outfn) // '.txt', action = 'write')
    write (upar, fmparec) "rho_0", a(c_rho,ngrid), "Central density"
    write (upar, fmparec) "temp_0", a(c_tgas,ngrid), "Central gas temperature"
    write (upar, fmparfc) "alpha", alpha, "Alpha parameter"
    write (upar, fmpare) "zscale", fzscale(mbh,mdot,radius)
    write (upar, fmpare) "radius", radius
    write (upar, fmparl) "has_corona", &
          & (cfg_temperature_method == EQUATION_BALANCE)
    write (upar, fmparl) "has_magnetic", .FALSE.
    write (upar, fmparl) "has_conduction", .FALSE.
    call wpar_gl(upar)
    close(upar)

    call write_results_3(z,ngrid,y,dy,y_labels,ny,a,a_labels,na,outfn)

    deallocate( z,y,dy,a )

contains

    subroutine rdconf(cfg)
        integer :: errno
        type(config), intent(inout) :: cfg
        character(len=2048) :: buf

        call mincf_get(cfg, "alpha", buf, errno)
        if ( iand(errno, mincf_not_found) .ne. 0 )  then
            error stop "Magnetic alpha-parameter (key: alpha) is REQUIRED!"
        end if
        read (buf,*) alpha

        call mincf_get(cfg, "radius", buf, errno)
        if ( iand(errno, mincf_not_found) .ne. 0 )  then
            error stop "Radius (key: radius) is REQUIRED!"
        end if
        read (buf,*) radius

    end subroutine

end program
