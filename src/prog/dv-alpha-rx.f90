program dv_alpha_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64
  use ieee_arithmetic, only: ieee_is_normal, ieee_is_nan
  use fileunits
  use relaxation
  use ss73solution, only: apxdisk, ss73_estim, ss73_refin
  use grid, only: space_linear, space_linlog

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model, errno
  integer :: ny = 3, i, globiter = 0
  integer, dimension(2) :: niter = [ 48, 12 ]
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:)
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, tau, y_trad
  real(dp) :: rhoc, Tc, Hdisk
  logical :: user_ff, user_bf
  integer, dimension(6) :: C_

  !----------------------------------------------------------------------------!

  ngrid = 576
  htop = 5

  !----------------------------------------------------------------------------!
  ! read command line

  parse_command_line: block

    integer :: i
    character(2**8) :: arg

    call rdargvgl

    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)

      ! use P_rad in alpha prescription?
      case ("-prad-alpha", "-alpha-prad")
        use_prad_in_alpha = .TRUE.
      case ("-no-prad-alpha", "-no-alpha-prad")
        use_prad_in_alpha = .FALSE.

      end select
    end do

  end block parse_command_line

  !----------------------------------------------------------------------------!
  ! read config

  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)

  user_ff = use_opacity_ff
  user_bf = use_opacity_bf

  !----------------------------------------------------------------------------!
  ! calculate the global parameters
  call cylinder(mbh, mdot, radius, rschw, omega, facc, teff, zscale)

  ! get the model number
  model = mrx_number( 'D', .FALSE., .FALSE. )
  call mrx_sel_nvar(model, ny)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid), M(ny*ngrid,ny*ngrid), tau(ngrid) )

  !----------------------------------------------------------------------------!

  call ss73_estim(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)
  call ss73_refin(mbh, mdot, radius, alpha, rhoc, Tc, Hdisk)

  !----------------------------------------------------------------------------!

  forall (i = 1:ngrid)
    x(i) = space_linlog(i,ngrid,htop) * Hdisk
    x0(i) = space_linlog(i,ngrid,htop) / htop
  end forall

  y_rho   => Y(C_(1)::ny)
  y_temp  => Y(C_(2)::ny)
  y_trad  => Y(C_(3)::ny)
  y_frad  => Y(C_(4)::ny)

  !----------------------------------------------------------------------------!

  forall (i = 1:ngrid)
    y_frad(i) = x0(i) * facc
    y_temp(i) = (1 - x0(i)) * (Tc - 0.841 * Teff) + 0.841 * Teff
    y_rho(i) =  rhoc * (exp(-0.5*(x(i)/Hdisk)**2) + 1e-8)
  end forall

  !----------------------------------------------------------------------------!

  use_opacity_ff = .false.
  use_opacity_bf = .false.

  call mrx_relax(model, x, y, niter(1), 0d0, errno)

  !----------------------------------------------------------------------------!

  use_opacity_ff = user_ff
  use_opacity_bf = user_bf

  if ( use_opacity_ff .or. use_opacity_bf ) then

    write (uerr,'(" --- ",A)') 'bf+ff opacity is ON'
    call mrx_relax(model, x, y, niter(2), 0.5d0, errno)

  end if

  !----------------------------------------------------------------------------!

  write (uerr,'(" --- ",A)') 'complete'

  open(33, file = trim(outfn) // '.dat', action = 'write')
  call wrdat(33)
  close(33)

  !----------------------------------------------------------------------------!


  open(33, file = trim(outfn) // '.txt', &
      action = 'write')
  call wrpar(33)
  call wpar_gl(33)
  close(33)

  open(34, file = trim(outfn) // ".col", action = 'write')
  call wrcol(34)
  close(34)

  deallocate(x,x0,Y,M,dY,tau)

contains

  subroutine wrdat(u)
    integer, intent(in) :: u
    integer :: i
    real(dp) :: pgas,prad
    call mktaues(x,y_rho,y_temp,tau)
    do i = 1,ngrid
      pgas =  cgs_k_over_mh / miu * y_rho(i) * y_temp(i)
      prad = cgs_a / 3 * y_trad(i)**4
      write (u,'(I7,*(ES14.5E3))') i, x(i), x(i) / zscale, tau(i), &
      y_rho(i), y_temp(i), y_frad(i),   &
      fksct(y_rho(i), y_temp(i)), fkabs(y_rho(i), y_temp(i)), &
      fkcnd(y_rho(i), y_temp(i)), pgas, prad
    end do
  end subroutine

  subroutine wrpar(u)
    integer, intent(in) :: u
    write (u, fmpari) "model", model
    write (u, fmpari) "niter", globiter
    write (u, fmpare) "radius", radius
    write (u, fmpare) "zscale", fzscale(mbh,mdot,radius)
    write (u, fmparec) "rho_0", rhoc, "Central density"
    write (u, fmparec) "temp_0", Tc, "Central gas temperature"
    write (u, fmparec) "Hdisk", Hdisk, "Disk height [cm]"
    write (u, fmparfc) "alpha", alpha, "Alpha parameter"
    write (u, fmparl) "has_corona", .FALSE.
    write (u, fmparl) "has_magnetic", .FALSE.
    write (u, fmparl) "has_conduction", .FALSE.
  end subroutine

  subroutine wrcol(u)
    integer, intent(in) :: u
    write(u, fmcol) 'i', 'i4'
    write(u, fmcol) 'z', 'f4'
    write(u, fmcol) 'h', 'f4'
    write(u, fmcol) 'tau', 'f4'
    write(u, fmcol) 'rho', 'f4'
    write(u, fmcol) 'temp', 'f4'
    write(u, fmcol) 'frad', 'f4'
    write(u, fmcol) 'ksct', 'f4'
    write(u, fmcol) 'kabs', 'f4'
    write(u, fmcol) 'kcnd', 'f4'
    write(u, fmcol) 'pgas', 'f4'
    write(u, fmcol) 'prad', 'f4'
  end subroutine

  !----------------------------------------------------------------------------!

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

  !----------------------------------------------------------------------------!

  subroutine mktaues(z, rho, T, tau)
    real(dp), intent(in), dimension(:) :: z, rho,T
    real(dp), intent(out), dimension(:) :: tau
    real(dp) :: dz, rhom, Tm
    integer :: i
    tau(size(tau)) = 1e-12
    do i = size(tau)-1, 1, -1
      rhom = (rho(i+1) + rho(i)) / 2
      Tm = (T(i+1) + T(i)) / 2
      dz = z(i) - z(i+1)
      tau(i) = tau(i+1) - kapes0(abuX,abuZ) * rhom * dz
    end do
  end subroutine

end program dv_alpha_relax
