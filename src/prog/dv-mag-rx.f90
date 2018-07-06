program dv_mag_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64
  use ieee_arithmetic, only: ieee_is_nan
  use fileunits
  use relaxation
  use slf_deriv, only: deriv
  use ss73solution, only: apxdisk, ss73_estim, ss73_refin
  use grid

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model, errno
  integer :: ny = 3, i, iter, nitert = 0, KL, KU
  integer, dimension(3) :: niter = [ 20, 8, 18 ]
  real(dp), allocatable, target :: x(:), x0(:), Y(:), dY(:), M(:,:), YY(:,:)
  real(dp), allocatable :: MB(:,:)
  real(dp), pointer :: yv(:,:)
  integer, dimension(:), allocatable :: ipiv
  real(dp) :: zphot, ztherm, zhard, zeqbc, ztmin, zcor
  logical, dimension(:), allocatable :: errmask
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, y_pmag, &
        & y_Trad, y_fcnd
  real(dp), allocatable, dimension(:) :: tau, heat
  real(dp) :: rho_0_ss73, temp_0_ss73, zdisk_ss73, err, err0, ramp, beta_0, qcor
  character(*), parameter :: fmiter = '(I5,2X,ES9.2,2X,F5.1,"%")'
  logical :: user_ff, user_bf, converged, has_corona
  integer, dimension(6) :: c_
  integer, parameter :: upar = 92
  !----------------------------------------------------------------------------!
  logical :: cfg_write_all_iters = .FALSE.
  character, parameter :: EQUATION_SIMPBALANCE = 'D'
  logical :: cfg_post_corona = .false.
  !----------------------------------------------------------------------------!

  real(dp), parameter :: typical_hdisk = 15

  integer, parameter :: ncols  = 34, &
      c_rho = 1, c_temp = 2, c_trad = 3, &
      c_pgas = 4, c_prad = 5, c_pmag = 6, &
      c_frad = 7, c_fmag = 8, c_fcnd = 9, &
      c_heat = 10, c_vrise = 11, &
      c_ksct = 12, c_kabs = 13, c_kabp = 14, &
      c_tau = 15, c_taues = 16, c_tauth = 17, &
      c_tavg = 18, c_beta = 19, c_coldens = 20, &
      c_kcnd = 21, c_coolb = 22, c_coolc = 23, &
      c_compy = 24, c_compy2 = 25, c_compfr = 26, c_fbfr = 27, &
      c_adiab1 = 28, c_gradad = 29, c_gradrd = 30, c_betamri = 31, &
      c_instabil = 32, c_qcor = 33, c_ionxi = 34

  character(8), dimension(ncols) :: labels

  labels(c_rho) = 'rho'
  labels(c_temp) = 'temp'
  labels(c_trad) = 'trad'
  labels(c_pgas) = 'pgas'
  labels(c_prad) = 'prad'
  labels(c_pmag) = 'pmag'
  labels(c_frad) = 'frad'
  labels(c_fmag) = 'fmag'
  labels(c_fcnd) = 'fcnd'
  labels(c_vrise) = 'vrise'
  labels(c_heat) = 'heat'
  labels(c_vrise) = 'vrise'
  labels(c_ksct) = 'ksct'
  labels(c_kabs) = 'kabs'
  labels(c_kabp) = 'kabp'
  labels(c_kcnd) = 'kcnd'
  labels(c_tau) = 'tau'
  labels(c_taues) = 'taues'
  labels(c_tauth) = 'tauth'
  labels(c_tavg) = 'tavg'
  labels(c_beta) = 'beta'
  labels(c_ionxi) = 'ionxi'
  labels(c_betamri) = 'betamri'
  labels(c_coldens) = 'coldens'
  labels(c_coolb) = 'coolb'
  labels(c_coolc) = 'coolc'
  labels(c_compy) = 'compy'
  labels(c_compy2) = 'compy2'
  labels(c_compfr) = 'compfr'
  labels(c_fbfr) = 'fbfr'
  labels(c_instabil) = 'instabil'
  labels(c_adiab1) = 'adiab1'
  labels(c_gradad) = 'gradad'
  labels(c_gradrd) = 'gradrd'
  labels(c_qcor) = 'qcor'

  !----------------------------------------------------------------------------!
  ! default values

  ngrid = -1
  htop = 180

  !----------------------------------------------------------------------------!
  ! initialize the globals, read the config etc

  call rdargvgl
  call rdargvrx
  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)


  open(upar, file = trim(outfn) // '.txt', action = 'write')

  user_ff = use_opacity_ff
  user_bf = use_opacity_bf
  converged = .true.
  has_corona = (cfg_temperature_method .ne. EQUATION_DIFFUSION) &
        .or. cfg_post_corona

  !----------------------------------------------------------------------------!
  ! check the magnetic parameters

  if (zeta <= 0 .or. zeta > 1) &
    error stop "zeta must be positive and less than 1"
  if (alpha <= 0) error stop "alpha must be positive"

  beta_0 = 2 * zeta / alpha + nu - 1
  qcor = 2 + alpha * (nu - 1) / zeta

  if (beta_0 < 0) error stop "MAGNETIC BETA cannot be negative! " &
        & // "zeta > alpha / 2 !!!"

  if (qcor < 1) error stop "I refuse to compute the disk where qcor < 1 !!!!"

  !----------------------------------------------------------------------------!
  ! some initial computation

  ! calculate the global parameters
  call cylinder(mbh, mdot, radius, rschw, omega, facc, teff, zscale)

  ! calculate the SS73 solution to obtain the initial profile for iteration
  call ss73_estim(mbh, mdot, radius, alpha, rho_0_ss73, temp_0_ss73, zdisk_ss73)
  call ss73_refin(mbh, mdot, radius, alpha, rho_0_ss73, temp_0_ss73, zdisk_ss73)

  !----------------------------------------------------------------------------!
  ! estimate the interval height

  if (cfg_auto_htop) then
    ! estimate the disk top from approximate solution from B15
    htop = (zdisk_ss73 / zscale) * sqrt((4 + alpha * nu / zeta) &
          * (1d-6**(-2 / (qcor + 2)) - 1))
    ! keep the disk dimension between 12H and 900H
    htop = min(max(htop, 12.0_dp), 900.0_dp)
  end if

  ! old method:
  ! htop = htop * (1 + 1 / beta_0) * max(0.4 * zdisk_ss73 / zscale, 1.0_dp)

  !----------------------------------------------------------------------------!
  ! if the grid number has not been set, choose the default

  if (ngrid .eq. -1) then
    select case (tgrid)
    case (grid_linear)
      ngrid = ceiling(90 * sqrt(htop))
    case (grid_log, grid_asinh)
      ngrid = ceiling(320 * log(1 + htop / typical_hdisk))
    case (grid_pow2)
      ngrid = ceiling(70 * sqrt(htop))
    case default
      error stop "this grid is not supported"
    end select
  end if

  !----------------------------------------------------------------------------!

  ! get the model number
  model = mrx_number( 'D', .TRUE., .FALSE. )
  call mrx_sel_nvar(model, ny)
  call mrx_sel_hash(model, C_)

  allocate( x(ngrid), x0(ngrid), Y(ny*ngrid), dY(ny*ngrid),  M(ny*ngrid,ny*ngrid), tau(ngrid), heat(ngrid), YY(ncols,ngrid) )
  allocate(errmask(ny*ngrid), ipiv(ny*ngrid))

  !----------------------------------------------------------------------------!
  ! generate the grid

  generate_grid: block

    select case (tgrid)
    case (grid_linear)
      forall (i = 1:ngrid)
        x(i)  = space_linear(i, ngrid, htop) * zscale
        x0(i) = space_linear(i, ngrid, htop) / htop
      end forall
    case (grid_log)
      forall (i = 1:ngrid)
        x(i)  = space_linlog(i, ngrid, htop / typical_hdisk) &
              * typical_hdisk * zscale
        x0(i) = space_linlog(i, ngrid, htop / typical_hdisk) &
              * typical_hdisk / htop
      end forall
    case (grid_asinh)
      forall (i = 1:ngrid)
        x(i)  = space_asinh(i, ngrid, htop / typical_hdisk) &
              * typical_hdisk * zscale
        x0(i) = space_asinh(i, ngrid, htop / typical_hdisk) &
              * typical_hdisk / htop
      end forall
    case (grid_pow2)
      forall (i = 1:ngrid)
        x(i)  = space_pow2(i, ngrid, real(100, dp)) * htop * zscale
        x0(i) = space_pow2(i, ngrid, real(100, dp))
      end forall
    case default
      error stop "this grid is not supported"
    end select

  end block generate_grid

  !----------------------------------------------------------------------------!
  ! set pointers

  y_rho   => Y(C_(1)::ny)
  y_temp  => Y(C_(2)::ny)
  y_trad  => Y(C_(3)::ny)
  y_frad  => Y(C_(4)::ny)
  y_pmag  => Y(C_(5)::ny)
  YV(1:ny,1:ngrid) => Y

  !----------------------------------------------------------------------------!
  ! generate the initial disk profile

  forall (i = 1:ngrid)
    y_frad(i) = x0(i) * facc
    y_temp(i) = (1 - x0(i)) * (temp_0_ss73 - 0.841 * Teff) + 0.841 * Teff
    y_rho(i) =  rho_0_ss73 * (exp(-0.5*(x(i)/zdisk_ss73)**2) + 1e-6)

    y_pmag(i) = 2 * cgs_k_over_mh * y_rho(i) * y_temp(i)   &
          & / (beta_0 * exp(- 0.5 * (x(i) / zdisk_ss73)**2 ) + 1e-2)
  end forall

  !----------------------------------------------------------------------------!
  ! do the initial relaxation with only electron scattering opacity

  if (cfg_write_all_iters) call saveiter(0)

  use_opacity_ff = .false.
  use_opacity_bf = .false.
  err0 = 0

  relx_opacity_es : do iter = 1,niter(1)

    call mrx_matrix(model, x, Y, M, dY)
    ! call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)
    call m2band(M, MB, KL, KU)
    call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

    forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
        & .and. ieee_is_normal(dY(i))
    err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
    ramp = ramp5(iter, niter(1) - 6)
    write(uerr,fmiter) nitert+1, err, 100*ramp

    Y(:) = Y + dY * ramp

    nitert = nitert + 1
    if (cfg_write_all_iters) call saveiter(nitert)

    if (err < 1e-5 .and. err0 / err > 5) then
      write (uerr, '("convergence reached with error = ",ES9.2)') err
      exit relx_opacity_es
    end if

    err0 = err

  end do relx_opacity_es


  !----------------------------------------------------------------------------!
  ! do relaxation with other opacities

  use_opacity_ff = user_ff
  use_opacity_bf = user_bf

  err0 = 0

  if ( use_opacity_bf .or. use_opacity_ff ) then

    write (uerr,*) '--- ff+bf opacity is on'

    relx_opacity_full : do iter = 1,niter(2)

      call mrx_matrix(model, x, Y, M, dY)
      ! call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)
      call m2band(M, MB, KL, KU)
      call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

      forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
          & .and. ieee_is_normal(dY(i))
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
      ramp = ramp2(iter, niter(2) - 2)
      write(uerr,fmiter) nitert+1, err, 100*ramp

      Y(:) = Y + dY * ramp

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

      if (err < 1e-5 .and. err0 / err > 5) then
        write (uerr, '("convergence reached with error = ",ES9.2)') err
        exit relx_opacity_full
      end if

      err0 = err

    end do relx_opacity_full


  end if

  !----------------------------------------------------------------------------!
  ! if the coorna was requested, relax the gas temperature

  if ( cfg_temperature_method .ne. EQUATION_DIFFUSION ) then

    write (uerr,*) '--- corona is on'

    err0 = 0

    call mrx_transfer(model, &
        mrx_number(cfg_temperature_method, .TRUE., use_conduction), x, Y)

    call mrx_sel_nvar(model, ny)
    call mrx_sel_hash(model, c_)

    deallocate(dY, M, errmask, ipiv)
    allocate(dY(ny*ngrid), M(ny*ngrid,ny*ngrid))
    allocate(errmask(ny*ngrid), ipiv(ny*ngrid))

    y_rho   => Y(C_(1)::ny)
    y_temp  => Y(C_(2)::ny)
    y_trad  => Y(C_(3)::ny)
    y_frad  => Y(C_(4)::ny)
    y_pmag  => Y(C_(5)::ny)
    YV(1:ny,1:ngrid) => Y

    if (cfg_temperature_method == EQUATION_BALANCE) niter(3) = niter(3) + 8

    if (use_conduction) then
      niter(3) = niter(3) + 12
      y_fcnd => Y(C_(6)::ny)
      y_fcnd(:) = 0
    end if

    call deriv(x, y_frad, heat)

    forall (i = 1:ngrid)
      y_temp(i) = heat(i) * (cgs_mel * cgs_c**2) &
            & / (16 * cgs_boltz * (cgs_kapes*y_rho(i)) &
            & * (cgs_stef*y_trad(i)**4) )
      y_temp(i) = sqrt(y_trad(i)**2 + y_temp(i)**2)
    end forall

    nitert = nitert + 1
    if (cfg_write_all_iters) call saveiter(nitert)

    !--------------------------------------------------------------------------!

    relx_corona : do iter = 1,niter(3)

      call mrx_matrix(model, x, Y, M, dY)
      ! call dgesv(size(M,2), 1, M, size(M,1), ipiv, dY, size(dY), errno)
      call m2band(M, MB, KL, KU)
      call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

      if (errno .ne. 0) exit relx_corona

      forall (i = 1:ngrid*ny) errmask(i) = (Y(i) .ne. 0) &
          & .and. ieee_is_normal(dY(i))
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
      ramp = ramp5(iter, niter(3) - 4) ! * min(err0 / err, 1.0_r64)

      write(uerr,fmiter) nitert+1, err, 100*ramp

      if (ieee_is_nan(err)) exit relx_corona

      if ((iter > 1 .and. err > err0 * 1.05) .or. (err > 1e4)) then
        write (uerr, '(''diverged: '', Es9.2, '' -> '', Es9.2)') err0, err
        converged = .false.
        exit relx_corona
      end if
      Y(:) = Y + dY * ramp

      where (y_trad < teff / 2) y_trad = teff / 2
      where (y_temp < teff / 2) y_temp = teff / 2
      where (.not. ieee_is_normal(y_rho) .or. y_rho < epsilon(1.0_dp) * rho_0_ss73) &
            y_rho = epsilon(1.0_dp) * rho_0_ss73

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

      if (err < 1e-5 .and. err0 / err > 5) then
        write (uerr, '("convergence reached with error = ",ES9.2)') err
        exit relx_corona
      end if

      err0 = err

    end do relx_corona

  end if

  write (uerr,*) '--- DONE'

  !----------------------------------------------------------------------------!
  ! write model data

  write_results: block

    open(33, file = trim(outfn) // '.dat', action = 'write')

    call fillcols(yv, c_, yy)

    do i = 1,ngrid
      write (33,'(I6,*(ES14.5E3))') i, x(i), x(i) / zscale, yy(:,i)
    end do

    close(33)

  end block write_results

  !----------------------------------------------------------------------------!
  ! write some global information

  call wpar_gl(upar)

  write (upar, fmhdr)  "disk information"
  write (upar, fmparfc) "alpha", alpha, "alpha parameter"
  write (upar, fmparfc) "eta", zeta, "field rise parameter"
  write (upar, fmparfc) "zeta", zeta, "field rise parameter"
  write (upar, fmparfc) "nu", nu, "reconnection parameter"
  write (upar, fmparfc) "qcor", qcor, "corona parameter"
  write (upar, fmparfc) "qcor_fact", maxval(yy(c_qcor,:)), "same incl. pressure"
  write (upar, fmparfc) "xcor", qcor, "corona parameter (old name)"

  write (upar, fmpare) "radius", radius
  write (upar, fmpare) "zscale", zscale
  write (upar, fmpare) "facc", facc
  write (upar, fmpare) "omega", omega

  associate (teff_real => (yy(c_frad, ngrid) / cgs_stef)**(0.25_dp))
    write (upar, fmpare) "teff", teff_real
    write (upar, fmparf) "teff_keV", teff_real * keV_in_kelvin
  end associate

  write (upar, fmpare) "rho_0", yy(c_rho,1)
  write (upar, fmpare) "temp_0", yy(c_temp,1)
  write (upar,fmpare) 'beta_0', yy(c_pgas,1) / yy(c_pmag,1)
  write (upar,fmpare) 'betakin_0', (yy(c_pgas,1) + yy(c_prad,1)) / yy(c_pmag,1)

  write (upar, fmpare) "rho_0_ss73", rho_0_ss73
  write (upar, fmpare) "temp_0_ss73", temp_0_ss73
  write (upar, fmparec) "zdisk_ss73", zdisk_ss73, "Disk height [cm]"
  write (upar, fmparfc) "hdisk_ss73", zdisk_ss73 / zscale, "Disk height [cm]"

  write (upar, fmhdr)  "Model information"
  write (upar, fmpari) "model", model
  write (upar, fmpari) "niter", nitert
  write (upar, fmpari) "ngrid", ngrid
  write (upar, fmparl) "converged", converged
  write (upar, fmparl) "has_corona", has_corona
  write (upar, fmparl) "has_magnetic", .TRUE.
  write (upar, fmparl) "has_conduction", use_conduction

  !----------------------------------------------------------------------------!
  ! save the column information for col2python

  write_columns: block

    open(35, file = trim(outfn) // ".col", action = 'write')

    write(35, fmcol) 'i', 'i4'
    write(35, fmcol) 'z', 'f4'
    write(35, fmcol) 'h', 'f4'

    do i = 1,ncols
      write (35,fmcol) labels(i), 'f4'
    end do

    close(35)

  end block write_columns

  !----------------------------------------------------------------------------!
  ! determine and save the photosphere location and other important values

  write_disk_globals: block

    use slf_interpol
    use slf_integrate
    real(dp) :: diskscale, tavgr
    logical :: has_eqbc = .false., has_tmin = .false.

    call save_interpolated(0.0_dp, 'midpl', 'midplane')

    call interpol(yy(c_tau,:), x, 0.667_dp, zphot)
    call save_interpolated(zphot, 'phot', 'tau = 1')

    call interpol(yy(c_tauth,:), x, 1.0_dp, ztherm)
    call save_interpolated(ztherm, 'therm', 'tau* = 1')

    if ( has_corona ) then
      zeqbc = -1
      call tabzero(x(ngrid:1:-1), yy(c_compfr,ngrid:1:-1), 0.5_dp, zeqbc)
      has_eqbc = zeqbc > 0
      write (upar, fmparl) 'has_eqbc', has_eqbc
      if (has_eqbc) then
        call save_interpolated(zeqbc, 'eqbc', 'compt. == brehms.')
      end if

      ztmin = -1
      call findtempmin(x, yy(c_temp,:), ztmin)
      has_tmin = ztmin > 0
      write (upar, fmparl) 'has_tmin', has_tmin
      if (has_tmin) then
        call save_interpolated(ztmin, 'tmin', 'temperature minimum')
      end if

      zcor = max(ztherm, ztmin)
      call save_interpolated(zcor, 'cor', 'max(ztherm, ztmin)')

    end if


    write (upar, fmhdr)  "some global parameters"
    write (upar, fmparfc) 'compy_0', yy(c_compy,1), 'total Y'
    write (upar, fmpare) 'frad_top', yy(c_frad, ngrid)
    write (upar, fmpare) 'fmag_top', yy(c_fmag, ngrid)
    write (upar, fmparf) 'fbfrac_top', yy(c_fbfr, ngrid)

    !--------------------------------------------------------------------------!
    ! instability location

    write(upar, fmparg) 'instabil', minval(yy(c_instabil,:))

    !--------------------------------------------------------------------------!
    ! compute the vertical disk scale and save it

    write (upar, fmhdr)  "column density and disk vertical scale"
    write (upar, fmparec) 'coldens', yy(c_coldens,ngrid), 'column density'

    ! vertical scale - weighted by density
    diskscale = sqrt(integrate(yy(c_rho,:) * x**2, x) / yy(c_coldens,ngrid))
    write (upar, fmparec) 'zdisk', diskscale, &
    'disk vertical scale (second moment)'
    write (upar, fmparf) 'hdisk', diskscale / zscale

    ! schwarzchild radius and disk ratio d = H / R
    write (upar, fmparg) 'ddisk', diskscale / (radius * rschw)

    ! average disk temperature
    tavgr = integrate(yy(c_rho,:) * yy(c_temp,:), x) / yy(c_coldens,ngrid)
    write (upar, fmparec) 'tavgr', tavgr, 'average temperature (by mass)'

    ! vertical scale - weighted by gas pressure
    diskscale = sqrt(integrate(yy(c_pgas,:) * x**2, x) &
                  /  integrate(yy(c_pgas,:),        x))
    write (upar, fmparec) 'zdisk_pgas', diskscale, &
         'disk height weighted by gas pressure'
    write (upar, fmparf) 'hdisk_pgas', diskscale / zscale

    ! vertical scale - weighted by magnetic pressure
    diskscale = sqrt(integrate(yy(c_pmag,:) * x**2, x) &
                  /  integrate(yy(c_pmag,:),        x))
    write (upar, fmparec) 'zdisk_pmag', diskscale, &
         'disk height weighted by magnetic pressure'
    write (upar, fmparf) 'hdisk_pmag', diskscale / zscale

    ! vertical scale - weighted by radiation pressure
    diskscale = sqrt(integrate(yy(c_prad,:) * x**2, x) &
                  /  integrate(yy(c_prad,:),        x))
    write (upar, fmparec) 'zdisk_prad', diskscale, &
         'disk height weighted by radiation pressure'
    write (upar, fmparf) 'hdisk_prad', diskscale / zscale

    write (upar, fmparf) 'hheat', integrate(yy(c_heat,:) * x, x) &
        / integrate(yy(c_heat,:), x) / zscale

  end block write_disk_globals

  !----------------------------------------------------------------------------!
  ! clean up

  close(upar)

  deallocate(x,x0,Y,M,dY,tau,heat,errmask,ipiv,yy)
  if (allocated(MB)) deallocate(MB)

  !----------------------------------------------------------------------------!

contains

  !----------------------------------------------------------------------------!
  ! save some important parameters, interpolated at given height.
  ! the results go to the .txt file
  ! usually, used to save properties at photosphere, base of the corona etc.

  subroutine save_interpolated(z, keyword, comment)

    use slf_interpol
    real(dp), intent(in) :: z
    character(*) :: keyword, comment
    real(dp) :: yz

    write (upar, fmhdr) 'properties in ' // comment

    write (upar, fmparec) 'z' // keyword, z, comment
    write (upar, fmparf) 'h' // keyword, z / zscale

    call interpol(x, yy(c_temp,:), z, yz)
    write (upar,fmparec) 'temp_' // keyword, yz, 'temperature in ' // comment
    write (upar,fmparf) 'temp_' // keyword // '_keV', yz * keV_in_kelvin

    call interpol(x, yy(c_rho,:), z, yz)
    write (upar,fmparec) 'rho_' // keyword, yz, 'density in ' // comment
    write (upar,fmparec) 'rho_' // keyword // '_nh', yz / cgs_mhydr, &
      'density in ' // trim(comment) // ' (nH)'

    call interpol(x, yy(c_trad,:), z, yz)
    write (upar,fmparec) 'trad_' // keyword, yz, 'rad temp in ' // comment
    write (upar,fmparf) 'trad_' // keyword // '_keV', yz * keV_in_kelvin

    ! average temperature of the corona
    yz = interpolf(x, yy(c_tavg,:), z) / interpolf(x, yy(c_tau,:), z)
    write (upar,fmparec) 'tavg_' // keyword, yz, 'average temperature in ' // comment
    write (upar,fmparf) 'tavg_' // keyword // '_keV', yz * keV_in_kelvin

    ! optical depth of the corona
    write (upar,fmparg) "tau_" // keyword, interpolf(x, yy(c_tau,:), z)
    write (upar,fmparg) "taues_" // keyword, interpolf(x, yy(c_taues,:), z)
    write (upar,fmparg) "tauth_" // keyword, interpolf(x, yy(c_tauth,:), z)

    write (upar,fmparf) 'compy_' // keyword, interpolf(x, yy(c_compy,:), z)
    write (upar,fmparf) 'compy2_' // keyword, interpolf(x, yy(c_compy2,:), z)

    write (upar,fmparec) 'compfr_' // keyword, interpolf(x, yy(c_compfr,:), z), &
          'compton / brehms. ratio'

    ! energy released in the corona
    associate (frad => interpolf(x, y_frad, z),       &
    &          frad_top => y_frad(ngrid),             &
    &          fmag => interpolf(x, yy(c_fmag,:), z))
      write (upar,fmpare) "frad_" // keyword, frad
      write (upar,fmparfc) "chi_" // keyword, 1 - frad / frad_top, &
      &   "relative amount of radiative flux released in the " // comment
      write (upar, fmparf) 'fbfrac_' // keyword, fmag / (frad + fmag)
    end associate

    ! magnetic beta
    write (upar,fmpargc) 'beta_' // keyword, interpolf(x, yy(c_beta,:), z), &
        & 'magnetic beta in ' // comment
    write (upar,fmparg) 'ionxi_' // keyword, interpolf(x, yy(c_ionxi,:), z)

    ! magnetic gradient
    write (upar,fmparg) 'qcor_' // keyword, interpolf(x, yy(c_qcor,:), z)

    ! column density: disk vs corona
    call interpol(x, yy(c_coldens,:), z, yz)
    associate (coldens => yy(c_coldens,ngrid))
      write (upar, fmparec) "coldens_below_" // keyword, yz, &
        "column density (disk only)"
      write (upar, fmparec) "coldens_" // keyword, coldens - yz, &
        "column density (corona only)"
      write (upar, fmparec) "mfrac_" // keyword, &
        (coldens - yz) / coldens, "mass fraction in " // comment
    end associate

    write (upar, fmparg) 'kapsct_' // keyword, interpolf(x, yy(c_ksct,:), z)
    write (upar, fmparg) 'kapabs_' // keyword, interpolf(x, yy(c_kabs,:), z)
    write (upar, fmparg) 'kapabp_' // keyword, interpolf(x, yy(c_kabp,:), z)
    write (upar, fmpargc) 'mag_gauss_' // keyword, &
    &  sqrt(8 * pi * interpolf(x, yy(c_pmag,:), z)), 'field strength in Gauss'

  end subroutine

  !----------------------------------------------------------------------------!
  ! save each relaxation iteration to a numbered file (001, 002, 003 etc)

  subroutine saveiter(iter)

    integer, intent(in) :: iter
    character(256) :: fn
    integer :: i

    write (fn,'(A,".",I0.3,".dat")') trim(outfn),iter
    open(33, file = trim(fn), action = 'write')

    call fillcols(yv, c_, yy)

    writeresults : do i = 1,ngrid
      write (33,'(I6,*(ES14.5E3))') i, x(i), x(i) / zscale, yy(:,i)
    end do writeresults

    close(33)

  end subroutine

  !----------------------------------------------------------------------------!

  subroutine fillcols(yv,c_,yy)
    procedure(funout_t), pointer :: fout
    real(dp) :: kabs,ksct,kabp,rhom,tempm,tradm,tcorrm,dx
    real(dp), dimension(:,:), intent(in) :: yv
    integer, dimension(:), intent(in) :: c_
    real(dp), dimension(:,:), intent(inout) :: yy
    integer :: i

    ! select the output function appropriate for this model
    call mrx_sel_fout(model, fout)

    ! evaluate the function for each point
    do i = 1,ngrid
      call fout(x(i), yv(:,i), yy(:,i))
    end do

    ! solve the exact balance after relaxation
    ! warning: this breaks strict hydrostatic equilibrium (but not much)
    if ( cfg_temperature_method /= EQUATION_BALANCE &
            .and. cfg_post_corona ) then
      post_corona: block
        use heatbalance, only: heatbil2
        real(dp) :: temp_old
        integer :: i
        do i = 1,ngrid
          temp_old = yy(c_temp,i)
          call heatbil2(yy(c_rho,i), yy(c_temp,i), yy(c_trad,i), &
                yy(c_heat,i), .false.)
          yy(c_pgas,i) = yy(c_pgas,i) * yy(c_temp,i) / temp_old
        end do
      end block post_corona
    end if

    ! opacities
    yy(c_ksct,:) = fksct(yy(c_rho,:), yy(c_temp,:))
    yy(c_kabs,:) = fkabs(yy(c_rho,:), yy(c_temp,:))
    yy(c_kabp,:) = fkabp(yy(c_rho,:), yy(c_temp,:))
    yy(c_kcnd,:) = fkcnd(yy(c_rho,:), yy(c_temp,:))

    ! cooling components: brehmstrahlung and compton
    cooling: block
      real(dp), dimension(ngrid) :: cb, cc, tcorr
      yy(c_coolb,:) = 4 * cgs_stef * yy(c_rho,:) * yy(c_kabp,:)   &
            * (yy(c_temp,:)**4 - yy(c_trad,:)**4)
      tcorr(:) = sqrt(1 + (4 * cgs_k_over_mec2 * yy(c_temp,:))**2)
      yy(c_coolc,:) = 4 * cgs_stef * yy(c_rho,:) * yy(c_ksct,:)   &
          * yy(c_trad,:)**4 * cgs_k_over_mec2 * 4 * (yy(c_temp,:) &
          * merge(tcorr(:), 1.0_dp, use_precise_balance) - yy(c_trad,:))
      cb(:) = yy(c_kabp,:) * yy(c_temp,:)**4
      cc(:) = yy(c_ksct,:) * yy(c_trad,:)**4 * cgs_k_over_mec2 &
          * 4 * yy(c_temp,:) * merge(tcorr(:), 1.0_dp, use_precise_balance)
      yy(c_compfr,:) = cc / (cb + cc)
    end block cooling

    ! radiative / total flux fraction
    yy(c_fbfr,1) = 0
    yy(c_fbfr,2:) = yy(c_frad,2:) / (yy(c_frad,2:) + yy(c_fmag,2:))

    ! radiation pressure fraction
    yy(c_ionxi,:) = yy(c_prad,:) / yy(c_pgas,:)

    ! integrate the optical depths and averaged temperature
    yy(c_tau,ngrid) = 0
    yy(c_taues,ngrid) = 0
    yy(c_tauth,ngrid) = 0
    yy(c_tavg,ngrid) = 0
    yy(c_compy,ngrid) = 0
    yy(c_compy2,ngrid) = 0

    integrate_tau: do i = ngrid-1,1,-1
      rhom = (yy(c_rho,i) + yy(c_rho,i+1)) / 2
      if (rhom < 0) rhom = 0
      tempm = (yy(c_temp,i) + yy(c_temp,i+1)) / 2
      tradm = (yy(c_trad,i) + yy(c_trad,i+1)) / 2
      dx = x(i+1) - x(i)

      kabs = merge(fkabs(rhom,tempm), 0.0_dp, tempm > 0)
      kabp = merge(fkabp(rhom,tempm), 0.0_dp, tempm > 0)
      ksct = merge(fksct(rhom,tempm), 0.0_dp, tempm > 0)

      yy(c_tau,  i) = yy(c_tau,  i+1) + dx * rhom * (kabs + ksct)
      yy(c_taues,i) = yy(c_taues,i+1) + dx * rhom * ksct
      yy(c_tauth,i) = yy(c_tauth,i+1) + dx * rhom * sqrt(kabp * (kabp + ksct))

      yy(c_tavg, i) = yy(c_tavg, i+1) + dx * rhom * (kabs + ksct) * tempm

      tcorrm = merge(sqrt(1 + (4 * cgs_k_over_mec2 * tempm)**2), 1.0_dp, &
        use_precise_balance)

      yy(c_compy, i) = yy(c_compy, i+1) + dx * rhom * ksct &
           * 4 * cgs_k_over_mec2 * (tempm * tcorrm - tradm)

      associate (tauesm => (yy(c_taues,i) + yy(c_taues,i+1)) / 2)
        yy(c_compy2, i) = yy(c_compy2, i+1) + dx * rhom * ksct &
          * 4 * cgs_k_over_mec2 * (tempm * tcorrm - tradm)  &
          * (2 * tauesm + 1)
          ! * (2 * tauesm**3 + tauesm) / sqrt(tauesm**4 + tauesm**2)
      end associate
    end do integrate_tau

    ! average tempearture
    ! yy(c_tavg,ngrid) = yy(c_temp,ngrid)
    ! yy(c_tavg,:ngrid-1) = yy(c_tavg,:ngrid-1) / yy(c_tau,:ngrid-1)

    ! magnetic beta parameter
    yy(c_beta,:) = yy(c_pgas,:) / yy(c_pmag,:)
    yy(c_betamri,:) = 2 * sqrt(yy(c_pgas,:) / yy(c_rho,:)) &
          / (omega * radius * rschw)

    ! here we compute d ln (cooling) / d ln T
    instability: block
      real(dp) :: coolcrit(ngrid)
      coolcrit(:) = 2 * cgs_stef * yy(c_rho,:) * yy(c_trad,:)**4  &
      * (   yy(c_kabp,:) * (9 - (yy(c_temp,:) / yy(c_trad,:))**4) &
      + 8 * yy(c_ksct,:) * cgs_k_over_mec2 * yy(c_temp,:))
      yy(c_instabil,:) = coolcrit(:) / yy(c_heat,:) - 1
    end block instability

    ! adiabatic gradients, according to Dalsgaard (book)
    gradients: block
      real(dp), dimension(ngrid) :: pgpt
      pgpt(:) = yy(c_pgas,:) / (yy(c_pgas,:) + yy(c_prad,:))
      yy(c_adiab1,:) = (32 - 24 * pgpt - 3 * pgpt**2) / (24 - 21 * pgpt)
      yy(c_gradad,:) = 2 * (4 - 3 * pgpt) / (32 - 24 * pgpt - 3 * pgpt**2)
      ! call deriv(log(yy(c_temp,:)), log(yy(c_pgas,:) + yy(c_prad,:)), &
            ! yy(c_gradrd,:))
      call loggrad(yy(c_temp,:), yy(c_pgas,:), yy(c_gradrd,:))
    end block gradients

    ! call loggrad(x, yy(c_pmag,:), yy(c_qcor,:))
    ! yy(c_qcor,:) = -yy(c_qcor,:)
    yy(c_qcor,:) = qcor - (yy(c_pgas,:) + yy(c_prad,:)) / yy(c_pmag,:) * (alpha / zeta)

    ! column density
    yy(c_coldens,1) = 0
    integrate_coldens: do i = 1, ngrid-1
      yy(c_coldens, i+1) = yy(c_coldens,i) &
      + (yy(c_rho,i) + yy(c_rho,i+1)) * (x(i+1) - x(i)) / 2
    end do integrate_coldens

  end subroutine

  !----------------------------------------------------------------------------!
  ! searches for temperature minimum in a given array

  subroutine findtempmin(x,temp,xtmin)
    real(dp), intent(in), dimension(:) :: x,temp
    real(dp), intent(inout) :: xtmin
    real(dp), dimension(size(x) - 1) :: dtemp,xm
    integer :: i

    do i = 1, size(x)-1
      xm(i) = (x(i) + x(i+1)) / 2
      dtemp(i) = temp(i+1) - temp(i)
    end do

    search_for_minimum: do i = 1, size(dtemp) - 1
      if (dtemp(i) .le. 0 .and. dtemp(i+1) .ge. 0) then
        xtmin = (dtemp(i+1)*xm(i) - dtemp(i)*xm(i+1)) / (dtemp(i+1) - dtemp(i))
        exit search_for_minimum
      end if
    end do search_for_minimum

  end subroutine

  !----------------------------------------------------------------------------!
  ! computes log gradient of any function

  pure subroutine loggrad(x, y, d)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: d(:)
    integer :: i

    if (size(x) /= size(y) .or. size(d) /= size(y)) error stop

    do concurrent (i = 2:size(y)-1)
      d(i) = (y(i+1) - y(i-1)) * (x(i+1) + x(i-1)) &
      &   / ((y(i+1) + y(i-1)) * (x(i+1) - x(i-1)) )
    end do

    d(1) = d(2)
    d(size(d)) = d(size(d) - 1)
  end subroutine

  !----------------------------------------------------------------------------!
  ! searches for zero in the array

  subroutine tabzero(x,y,y0,x0)
    real(dp), intent(in) :: x(:), y(:), y0
    real(dp), intent(inout) :: x0
    integer :: i

    if (size(x) /= size(y)) error stop "tabzero: size(x) /= size(y)"

    search_for_zero: do i = 1, size(y) - 1
      if ((y(i) - y0) * (y(i+1) - y0) .le. 0) then
        x0 = ((y(i+1) - y0) * x(i) - (y(i) - y0) * x(i+1)) / (y(i+1) - y(i))
        exit search_for_zero
      end if
    end do search_for_zero

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

    call mincf_get(cfg, "eta", buf, errno)
    if ( iand(errno, mincf_not_found) .ne. 0 )  then
      error stop "Field rise speed (key: eta) is REQUIRED!"
    end if
    read (buf,*) zeta

    call mincf_get(cfg, "nu", buf, errno)
    if ( iand(errno, mincf_not_found) .eq. 0 )  then
      read (buf,*) nu
    else
      write (0,*) "warning: assuming nu = 0"
      nu = 0
    end if

  end subroutine

  !----------------------------------------------------------------------------!

  subroutine rdargvrx
    integer :: i
    character(2**8) :: arg

    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)

      ! write every iteration to another file? useful to demonstrate relaxation
      case ("-write-all","-all")
        cfg_write_all_iters = .TRUE.

      ! recalculate the cooling-heating balance after relaxation? works best
      ! with -compton switch or alone (not much sense with -corona switch)
      case ("-post-corona")
        cfg_post_corona = .TRUE.
      case ("-no-post-corona")
        cfg_post_corona = .FALSE.

      ! include relativictic term in Compton source function? may cause
      ! some inconsistencies.
      case ("-relativistic", "-rel", "-relcompt")
        use_precise_balance = .TRUE.
      case ("-no-relativistic", "-no-rel", "-no-relcompt")
        use_precise_balance = .FALSE.

      ! enable PP condition for MRI shutdown? can cause trouble for convergence
      case ("-quench","-quench-mri","-qmri")
        use_quench_mri = .TRUE.
      case ("-no-quench","-no-quench-mri","-no-qmri")
        use_quench_mri = .FALSE.

      ! use P_rad in alpha prescription?
      case ("-prad-alpha", "-alpha-prad")
        use_prad_in_alpha = .TRUE.
      case ("-no-prad-alpha", "-no-alpha-prad")
        use_prad_in_alpha = .FALSE.

      end select
    end do
  end subroutine

  !----------------------------------------------------------------------------!

end program
