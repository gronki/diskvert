program dv_mag_relax

  use confort
  use globals
  use settings
  use iso_fortran_env, only: sp => real32, dp => real64, int64
  use ieee_arithmetic, only: ieee_is_nan, ieee_is_normal
  use fileunits
  use relaxation
  use slf_deriv, only: deriv
  use ss73solution, only: apxdisk, ss73_estim, ss73_refin
  use grid
  use slf_interpol

  !----------------------------------------------------------------------------!
  implicit none

  type(config) :: cfg
  integer :: model, errno
  integer :: ny = 3, nitert = 0
  real(dp), allocatable, target :: x(:),  Y(:), YY(:,:)
  real(dp), pointer :: yv(:,:), dyv(:,:)
  real(dp), pointer, dimension(:) :: y_rho, y_temp, y_frad, y_pmag, &
        & y_Trad, y_fcnd
  real(dp) :: rho_0_ss73, temp_0_ss73, zdisk_ss73, qcor, beta_0
  character(*), parameter :: fmiter = '(I5,2X,ES9.2,2X,F5.1,"%  ")'
  character(*), parameter :: fmiterw = '("' // achar(27) // '[33;1m",I5,2X,ES9.2,2X,F5.1,"%  ' // achar(27) // '[0m")'
  logical :: user_ff, user_bf, converged, has_corona
  integer, dimension(6) :: c_
  integer, parameter :: upar = 92
  !----------------------------------------------------------------------------!
  integer(int64) :: timing(2)
  logical :: with_perf = .false.
  !----------------------------------------------------------------------------!
  logical :: cfg_write_all_iters = .FALSE.
  character, parameter :: EQUATION_SIMPBALANCE = 'D'
  logical :: cfg_smooth_delta = .false.
  logical :: cfg_post_corona = .false., cfg_simple_post = .false., &
    cfg_magnetic = .true., cfg_trim_vacuum = .false.
  character(len=8) :: cfg_corestim_method = ""
  !----------------------------------------------------------------------------!
  real(dp), parameter :: trim_density_thresh = 1e-11
  real(dp), parameter :: typical_hdisk = 12
  real(dp) :: grid_nonlnr = 1e-3
  !----------------------------------------------------------------------------!


  integer, parameter :: ncols  = n_yout + 35, &
      c_ksct      = n_yout +  1, &
      c_kabs      = n_yout +  2, &
      c_kabp      = n_yout +  3, &
      c_tau       = n_yout +  4, &
      c_taues     = n_yout +  5, &
      c_tauth     = n_yout +  6, &
      c_tavg      = n_yout +  7, &
      c_beta      = n_yout +  8, &
      c_coldens   = n_yout +  9, &
      c_kcnd      = n_yout + 10, &
      c_coolb     = n_yout + 11, &
      c_coolc     = n_yout + 12, &
      c_compy     = n_yout + 13, &
      c_compy2    = n_yout + 14, &
      c_compfr    = n_yout + 15, &
      c_fbfr      = n_yout + 16, &
      c_adiab1    = n_yout + 17, &
      c_gradad    = n_yout + 18, &
      c_gradrd    = n_yout + 19, &
      c_pmagmri   = n_yout + 20, &
      c_instabil  = n_yout + 21, &
      c_qcor      = n_yout + 22, &
      c_ionxi     = n_yout + 23, &
      c_dnh       = n_yout + 24, &
      c_nhtot     = n_yout + 25, &
      c_betagen   = n_yout + 26, &
      c_cool      = n_yout + 27, &
      c_cool_dr   = n_yout + 28, &
      c_cool_dT   = n_yout + 29, &
      c_cool_ref  = n_yout + 30, &
      c_rho_max   = n_yout + 31, &
      c_radpz     = n_yout + 32, &
      c_cf_x      = n_yout + 33, &
      c_cf_y      = n_yout + 34, &
      c_instabv   = n_yout + 35

  !----------------------------------------------------------------------------!

  character(12), dimension(ncols) :: labels

  labels(c_rho) = 'rho'
  labels(c_temp) = 'temp'
  labels(c_trad) = 'trad'
  labels(c_pgas) = 'pgas'
  labels(c_prad) = 'prad'
  labels(c_pmag) = 'pmag'
  labels(c_frad) = 'frad'
  labels(c_fmag) = 'fmag'
  labels(c_fcnd) = 'fcnd'
  labels(c_ptot_gen) = 'ptot_gen'
  labels(c_heat) = 'heat'
  labels(c_heatm) = 'heatm'
  labels(c_heatr) = 'heatr'
  labels(c_vrise) = 'vrise'
  labels(c_qmri) = 'qmri'

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
  labels(c_pmagmri) = 'pmagmri'
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
  labels(c_dnh) = 'dnh'
  labels(c_nhtot) = 'nhtot'
  labels(c_betagen) = 'betagen'
  labels(c_cool) = 'cool'
  labels(c_cool_dr) = 'cool_dr'
  labels(c_cool_dT) = 'cool_dT'
  labels(c_cool_ref) = 'cool_ref'
  labels(c_rho_max) = 'rho_max'
  labels(c_radpz) = 'radpz'
  labels(c_cf_x) = 'cf_x'
  labels(c_cf_y) = 'cf_y'
  labels(c_instabv) = 'instabv'

  !----------------------------------------------------------------------------!

  ngrid = -1 

  !----------------------------------------------------------------------------!
  ! initialize the globals, read the config etc

  write(*, '("diskvert v.", a)') version

  call rdargvgl
  call rdargvrx
  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)

  write(*, '(a, 1x, a)') "FILE", trim(outfn)
  open(upar, file = trim(outfn) // '.txt', action = 'write')

  !----------------------------------------------------------------------------!

  user_ff = use_opacity_ff
  user_bf = use_opacity_bf
  converged = .true.
  has_corona = cfg_magnetic .and. &
    ((cfg_temperature_method .ne. EQUATION_DIFFUSION) .or. cfg_post_corona)

  !----------------------------------------------------------------------------!

  if (with_perf) call system_clock(timing(1))

  !----------------------------------------------------------------------------!
  ! some initial computation

  ! calculate the global parameters
  call cylinder(mbh, mdot, radius, rschw, omega, facc, teff, zscale)

  ! calculate the SS73 solution to obtain the initial profile for iteration
  call ss73_estim(mbh, mdot, radius, alpha, rho_0_ss73, temp_0_ss73, zdisk_ss73)
  call ss73_refin(mbh, mdot, radius, alpha, rho_0_ss73, temp_0_ss73, zdisk_ss73)

  !----------------------------------------------------------------------------!
  ! check the magnetic parameters

  if (alpha <= 0) error stop "alpha must be positive"

  if (cfg_magnetic) then
    if (eta > 1) write (0, *) "attention: eta > 1"
    if (eta <= 0) error stop "eta <= 0"

    beta_0 = (2 * eta + alpha * (nu - 1)) / alpha

    write(*, '(a10, f10.2)') 'beta_0 = ', beta_0
    if (beta_0 < 0) error stop "beta_0 < 0"

    qcor = (2 * eta + alpha * (nu - 1)) / eta
    write(*, '(a10, f10.2)') 'qcor = ', qcor
    write(*, '(a10, f10.2)') 'A-a = ', (2 * eta + alpha * (nu - 1))

    if (qcor <= 1) write (0, *) 'warning: qcor <= 1'
    ! if (qcor <= 1 .and. use_flux_correction .and. .not. use_quench_mri) &
    ! &   error stop "qcor <= 1 .and. use_flux_correction .and. .not. use_quench_mri"
    if (qcor <= 0) error stop "qcor <= 0"
  end if

  !----------------------------------------------------------------------------!
  ! estimate the interval height

  if (cfg_auto_htop) then
    block_auto_height: block
      real(dp) :: hg, hr, hm
      real(dp), parameter :: pmag_cutoff = 1e-8

      hg = zdisk_ss73 / zscale
      hr = cgs_kapes * facc / (omega**2 * cgs_c * zscale)
      write(*, '("SS73 height    ", f10.1)') hg
      write(*, '("rad press heigh", f10.1)') hr

      if (cfg_magnetic) then
        ! hm = sqrt((4 + alpha * nu / eta) * (pmag_cutoff**(-2 / (qcor + 1)) - 1))
        hm = 1 + pmag_cutoff**(-1 / (qcor + 2))
        write(*, '("magnetic height", f10.1)') hm
        htop = (hg + hr) * hm
        htop = min(max(htop, 9._dp), 3e4_dp)
      else
        htop = 9 * hg + hr
      end if
      
      write(*, '("assumed height ", f10.1)') htop
    end block block_auto_height
  end if

  !----------------------------------------------------------------------------!
  ! if the grid number has not been set, choose the default

  grids_block: block
    real(dp) :: a

    a = qcor / 20
    a = log(1 + qcor) / log(1 + 50._dp)
    a = max(0._dp, min(1._dp, a))

    if (ngrid == -1) then
      if (cfg_magnetic) then
        ngrid = 800 + nint((2500 - 800) * a)
      else
        ngrid = 1024
      end if
    end if

    ngrid = nint(ngrid / 32.) * 32
    if (cfg_magnetic) grid_nonlnr = 10**(-3.5 + 3.0 * a)

    write(*, '("GRID     a = ", f4.2)') a
    write(*, '("GRID ngrid = ", i4)') ngrid
    write(*, '("GRID nonln = ", f4.1)') log10(grid_nonlnr)
  end block grids_block


  !----------------------------------------------------------------------------!

  ! get the model number
  model = mrx_number( 'D', cfg_magnetic, use_conduction )
  call mrx_sel_nvar(model, ny)
  call mrx_sel_hash(model, C_)

  allocate(x(ngrid), Y(ny*ngrid), YY(ncols,ngrid))

  !----------------------------------------------------------------------------!
  ! generate the grid

  call generate_grid(x, htop * zscale)

  !----------------------------------------------------------------------------!
  ! set pointers

  YV(1:ny,1:ngrid) => Y
  y_rho   => yv(c_(1),:)
  y_temp  => yv(c_(2),:)
  y_trad  => yv(c_(3),:)
  y_frad  => yv(c_(4),:)
  if (cfg_magnetic) y_pmag => yv(C_(5),:)

  !----------------------------------------------------------------------------!
  ! generate the initial disk profile

  initial_profile: block
    integer :: i
    real(dp) :: pcentr, pth, x0(ngrid)

    x0(:) = x / x(ngrid)

    do i = 1, ngrid
      y_frad(i) = ramp6r(min(x0(i) / 0.2, 1.0_dp)) * facc
      y_temp(i) = (1 - x0(i))**2 * (temp_0_ss73 - 0.841 * Teff) + 0.841 * Teff
      y_rho(i) =  2 * rho_0_ss73 * exp(-0.5*(x(i) / zdisk_ss73)**2)

      if (cfg_magnetic) then
        y_rho(i) = y_rho(i) + rho_0_ss73 * 0.1 / (1 + (x(i) / zdisk_ss73)**(qcor + 2))
      else
        y_rho(i) = y_rho(i) + rho_0_ss73 * 1e-8
      end if

      if (cfg_magnetic) then
        pcentr = cgs_k_over_mh * temp_0_ss73 * rho_0_ss73 / miu &
        & + merge((cgs_a / 3) * temp_0_ss73**4, 0._dp, use_prad_in_alpha)
        y_pmag(i) = pcentr / beta_0 * (1 + (0.5 * x(i) / zdisk_ss73)**2)**(-qcor/2)
        pth = cgs_k_over_mh * y_temp(i) * y_rho(i) / miu &
        & + merge((cgs_a / 3) * y_temp(i)**4, 0._dp, use_prad_in_alpha)
        ! y_pmag(i) = max(y_pmag(i), 2 * pth / (beta_0 + nu))
      end if
    end do

    if (cfg_write_all_iters) call saveiter(0)

  end block initial_profile

  !----------------------------------------------------------------------------!

  relaxation_block: block

    real(dp), allocatable, target :: dY(:), M(:,:)
    real(dp), allocatable :: MB(:,:)
    integer, dimension(:), allocatable :: ipiv
    logical, dimension(:), allocatable :: errmask
    real(dp) :: err, err0, ramp, rholim, ers
    integer :: iter, kl, ku, i, negative_rho !, iter_opacity

    allocate(errmask(ny*ngrid), ipiv(ny*ngrid), dY(ny*ngrid), M(ny*ngrid,ny*ngrid))
    M(:,:) = 0
    dyv(1:ny, 1:ngrid) => dY

    !----------------------------------------------------------------------------!
    ! do the initial relaxation with only electron scattering opacity

    err0 = -1
    ! iter_opacity = 0
    ramp = 0
    qmri_kill = 0.0
    threshpow = 1.0
    opacities_kill = 0.
    negative_rho = 0

    converged = .false.

    relx_disk : do iter = 1, 2000

      if (any(y_rho <= 0)) then
        negative_rho = negative_rho + 1
      else
        negative_rho = 0
      end if

      if ((opacities_kill <= 0 .and. iter > 5 .and. err < 0.03 .and. .not. negative_rho > 0) &
      &   .or. (opacities_kill > 0 .and. err < 0.3)) then
        opacities_kill = min(1._dp, (opacities_kill + 0.25) * 1.05 - 0.25)
      end if

      if (iter > 500 .and. err < 1e-5 .and. negative_rho > 0) exit relx_disk

      call mrx_matrix(model, x, Y, M, dY)
      call mrx_bandim(model, kl, ku)
      call m2band(M, KL, KU, MB)
      call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

      errmask(:) = (Y .ne. 0) .and. ieee_is_normal(dY)
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))

      if (ieee_is_nan(err) .or. err > 1e12) then
        write(*, '("'// achar(27) //'[1;31;7mdiverged: ", Es9.2, " -> ", Es9.2, "'&
            // achar(27) //'[0m")') err0, err
        exit relaxation_block
      end if

      ers = merge((err / err0)**0.2, 1._dp, err0 > 0)
      ers = max(0.5_dp, min(4.0_dp, ers))
      ramp = max(err_ramp(err * ers, 0.5_dp, 0.6_dp), 1e-3_dp)

      ! rholim = minval(-yv(:,c_(1)) / dyv(:,c_(1)), &
      ! & mask=yv(:,c_(1)) > 0 .and. dyv(:,c_(1)) < 0 .and. ieee_is_normal(-yv(:,c_(1)) / dyv(:,c_(1))))
      ! rholim = min(1._dp, rholim)
      ! if (rholim < 1) ramp = ramp * min(1._dp, max(0.1_dp, 0.8 * rholim / ramp))
        
      write(*, '(a, i5, 1x, 1es10.2, 4x, "ramps=", 3(2x, F5.1, "%"), 4x, a, a)') &
      & trim(merge(achar(27) // '[33;1m', repeat(' ', 7), err0 < err .and. err0 > 0)), &
      & nitert+1, err,  100*ramp, 100*opacities_kill, 100*qmri_kill,  &
      & trim(merge('rho < 0', repeat(' ', 7), negative_rho > 0)), &
      & trim(merge(achar(27) // '[0m', repeat(' ', 4), err0 < err .and. err0 > 0))

      Y(:) = Y + dY * ramp

      ! fix negative densities
      if (iter > 20 .and. negative_rho > 5) then
        where (y_rho < 0) y_rho = -y_rho
        call smooth2(y_rho, 3)
        where (y_rho < 0) y_rho = -y_rho
        call smooth2(y_rho, 3)
        where (y_rho < 0) y_rho = -y_rho
        call smooth2(y_rho, 3)
      end if

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)

      if (err < 0.03 .and. use_quench_mri .and. .not. negative_rho > 0) then
        qmri_kill = min(1._dp, (qmri_kill + 0.25) * 1.03 - 0.25)
        threshpow =  1 + 3 * qmri_kill
      end if

      if (err < 1e-5 .and. (err0 / err > 5 .or. err < 1e-8) .and. opacities_kill > 0.999 &
      &   .and. (qmri_kill > 0.999 .or. .not. use_quench_mri) .and. negative_rho <= 0) then
        write(*, '("STRUCTURE convergence reached with error = '// achar(27) &
            // '[1m",ES9.2,"'// achar(27) //'[0m")') err
        converged = .true.
        exit relx_disk
      end if

      err0 = err

      if (cfg_trim_vacuum .and. (iter > 5 .and. opacities_kill > 0.3) &
      & .and. any(y_rho(ngrid/2:) < trim_density_thresh * maxval(y_rho))) then
        trim_space_vacuum: block
          use slf_interpol, only: interpol

          real(dp) :: xcut
          real(dp), allocatable :: xcopy(:), ycopy(:)
          integer :: i,j

          call tabzero(x, y_rho, 2 * trim_density_thresh * y_rho(1), xcut)
          if (xcut / zscale < 3 .or. xcut / x(ngrid) > 0.98) exit trim_space_vacuum

          write(*, '("'// achar(27) //'[96;1m<<< trimming top from ", &
          &   f6.1, " to ", f6.1, "'// achar(27) //'[0m")') &
          &   x(ngrid) / zscale, xcut / zscale

          xcopy = x(:)
          ycopy = y(:)

          call generate_grid(x, xcut)

          do j = 1, ny
            associate (y_j => y(j::ny), ycopy_j => ycopy(j::ny))
              do i = 1, ngrid
                call interpol(xcopy, ycopy_j, x(i), y_j(i))
              end do
            end associate
          end do
        end block trim_space_vacuum
      end if

    end do relx_disk

    if (opacities_kill < 0.99 .or. err > 5e-2) then
      write(*, '("'// achar(27) //'[1;31;7mnot converged'&
          // achar(27) //'[0m")')
      converged = .false.
      exit relaxation_block
    end if

    !----------------------------------------------------------------------------!
    ! if the coorna was requested, relax the gas temperature

    if (cfg_magnetic .and. cfg_temperature_method .ne. EQUATION_DIFFUSION) then

      write(*,*) '--- corona is on'

      err0 = 0

      call mrx_transfer(model, &
      mrx_number(cfg_temperature_method, cfg_magnetic, use_conduction), x, Y)

      call mrx_sel_nvar(model, ny)
      call mrx_sel_hash(model, c_)

      deallocate(dY, M, errmask, ipiv)
      allocate(dY(ny*ngrid), M(ny*ngrid,ny*ngrid))
      M(:,:) = 0
      allocate(errmask(ny*ngrid), ipiv(ny*ngrid))
      
      yv(1:ny, 1:ngrid) => Y
      dyv(1:ny, 1:ngrid) => dY
      y_rho  => yv(c_(1),:)
      y_temp => yv(c_(2),:)
      y_trad => yv(c_(3),:)
      y_frad => yv(c_(4),:)
      y_pmag => yv(c_(5),:)

      if (use_conduction) then
        y_fcnd => yv(c_(6),:)
        y_fcnd(:) = 0
      end if

      !--------------------------------------------------------------------------!


      estimate_corona: block
        use heatbalance, only: heatbil2
        real(dp) :: heat(ngrid), tmax, wsp1, tmax2
        real(dp), parameter :: mixin = 0.5
        integer :: i

        ! y_pgas(:) = cgs_k_over_mh / miu * y_rho(:) * y_trad(:)
        ! y_prad(:) = cgs_a * y_trad(:)**4 / 3
        ! heat(:) = 2 * (eta + alpha * nu) * omega * y_pmag(:) &
        !   - alpha * omega * (y_pgas(:) + y_pmag(:) &
        !   + merge(y_prad(:), 0.0_dp, use_prad_in_alpha))
        heat(:) = (2 * eta + alpha * (2 * nu - 1)) * omega * y_pmag(:)

        if (cfg_corestim_method == '') then
          select case(cfg_temperature_method)
          case (EQUATION_BALANCE)
            cfg_corestim_method = 'full'
          case (EQUATION_COMPTON)
            cfg_corestim_method = 'compt'
          end select
        end if

        select case(cfg_corestim_method)

        case ('compt')
          do concurrent (i = 1:ngrid)
            y_temp(i) = heat(i) * (cgs_mel * cgs_c**2) &
              & / (16 * cgs_boltz * (cgs_kapes * y_rho(i)) &
              & * (cgs_stef * y_trad(i)**4) )
            y_temp(i) = y_trad(i)**(1 - mixin) * (y_trad(i) + y_temp(i))**mixin
          end do

        case ('full')
          do concurrent (i = 1:ngrid)
            y_temp(i) = y_trad(i)
            call heatbil2(y_rho(i), y_temp(i), y_trad(i), heat(i), .false.)
            y_temp(i) = y_trad(i)**(1 - mixin) * y_temp(i)**mixin
          end do

        case ('instab')
          wsp1 = (8. / 3.) * cgs_kapes / kram0p(abuX, abuZ) * cgs_k_over_mec2
          do concurrent (i = 1:ngrid)
            y_temp(i) = y_trad(i)
            call heatbil2(y_rho(i), y_temp(i), y_trad(i), heat(i), .false.)
            tmax = 1.38 * y_trad(i) + (wsp1 * y_trad(i)**5 / y_rho(i))**2
            ! tmax2 = 1.2 * (2 * eta / alpha + 2 * nu - 1) * y_pmag(i) / (cgs_k_over_mh * y_rho(i) / 0.5)
            y_temp(i) = min(y_temp(i), tmax)
            ! y_temp(i) = 1.06 / (1 / y_temp(i)**4 + 1 / tmax**4)**(1. / 4.) 
            y_temp(i) = y_trad(i)**(1 - mixin) * y_temp(i)**mixin
          end do
          call smooth2(y_temp, 3)

        case default
          y_temp(:) = y_trad
          write(*, *) 'warning: no corona estimation'
        end select

      end block estimate_corona

      !--------------------------------------------------------------------------!

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)
      
      !--------------------------------------------------------------------------!
      
      converged = .false.

      relx_corona : do iter = 1, 300

        call mrx_matrix(model, x, Y, M, dY)
        call mrx_bandim(model, kl, ku)
        call m2band(M, KL, KU, MB)
        call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

        if (errno .ne. 0) then
          write(*, *) 'DGBSV ERRNO=', errno
          exit relx_corona
        end if

        if (cfg_smooth_delta) then
          ! where (.not. ieee_is_normal(dyv)) dyv = 0
          do concurrent (i = 1:ny)
            call smooth2(dYv(i,:), 2)
          end do
        end if
        
        errmask(:) = (Y .ne. 0) .and. ieee_is_normal(dY)
        err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
        ramp = max(err_ramp(err, 0.10_dp, 0.50_dp), 1e-4_dp)
        
        if (iter > 1 .and. err > err0) then
          write(*,fmiterw) nitert+1, err, 100*ramp
        else
          write(*,fmiter) nitert+1, err, 100*ramp
        end if
        
        where (.not. ieee_is_normal(dyv)) dyv = 0
        Y(:) = Y + dY * ramp
        where (.not. ieee_is_normal(Y)) Y = 0
        
        ! call smooth2(y_temp, 3)
        ! call smooth2(y_rho, 2)

        nitert = nitert + 1
        if (cfg_write_all_iters) call saveiter(nitert)

        if (ieee_is_nan(err) .or. (iter > 1 .and. err > err0 * 300) .or. (err > 1e9)) then
          write(*, '("' // achar(27) // '[1;31;7mdiverged: ", Es9.2, " -> ", Es9.2, "'&
              // achar(27) //'[0m")') err0, err
          if (iter > 2) write(*, '(a, i4)') 'ncorit = ', iter
          exit relaxation_block
        end if

        if (err < 5e-5 .and. (err0 / err > 8 .or. err < 1e-8)) then
          write(*, '("CORONA convergence reached with error = '// achar(27) &
              // '[1m",ES9.2,"'// achar(27) //'[0m")') err
          converged = .true.
          exit relx_corona
        end if

        err0 = err
      end do relx_corona
    end if

    ! solve the exact balance after relaxation
    ! warning: this breaks strict hydrostatic equilibrium (but not much)
    if (converged .and. cfg_post_corona .and. (.not. cfg_simple_post)) then
      put_post_corona: block

      use heatbalance, only: heatbil2
      integer :: i
      real(dp) :: temp_old

      if (cfg_temperature_method == EQUATION_DIFFUSION) then
        call mrx_transfer(model, &
          mrx_number(EQUATION_BALANCE, cfg_magnetic, use_conduction), x, Y)

        call mrx_sel_nvar(model, ny)
        call mrx_sel_hash(model, c_)

        deallocate(dY, M, errmask, ipiv)
        allocate(dY(ny*ngrid), M(ny*ngrid,ny*ngrid))
        M(:,:) = 0
        allocate(errmask(ny*ngrid), ipiv(ny*ngrid))
        
        yv(1:ny, 1:ngrid) => Y
        dyv(1:ny, 1:ngrid) => dY
        y_rho  => yv(c_(1),:)
        y_temp => yv(c_(2),:)
        y_trad => yv(c_(3),:)
        y_frad => yv(c_(4),:)
        y_pmag => yv(c_(5),:)

        if (use_conduction) then
          y_fcnd => yv(c_(6),:)
          y_fcnd(:) = 0
        end if
      end if

      ! evaluate model because we are going to require heating
      call evaluate_fout(model, x, yv, yy)

      do iter = 1, 33

        do concurrent (i = 1:ngrid)
          temp_old = yy(c_temp,i)
          call heatbil2(yy(c_rho,i), yy(c_temp,i), yy(c_trad,i), &
                yy(c_heat,i), .false.)
          yv(c_(2),i) = yy(c_temp,i)**0.4 * temp_old**0.6
        end do

        err = sum((2 * (yv(c_(2),:) - yy(c_temp,:)) / (yv(c_(2),:) + yy(c_temp,:)))**2)
        write(*, '(a, 1x, i3, 1x, es8.1)') 'POST', iter, err
    
        call evaluate_fout(model, x, yv, yy)

        nitert = nitert + 1
        if (cfg_write_all_iters) call saveiter(nitert)

        if (any(ieee_is_nan(yv))) then
          write(*, '(a, i0)') 'post corona fail at iter=', iter
          converged = .false.
          exit put_post_corona
        else if (err < 1e-2) then
          exit put_post_corona
        end if
        
      end do

      write(*, '(a)') 'POST CORONA OK'

      end block put_post_corona
    end if 

    if (converged) write(*, '(2X,a)') achar(27) // '[1;7;32m   SUCCESS   ' // achar(27) // '[0m'

  end block relaxation_block


  !----------------------------------------------------------------------------!
  ! write model data

  call saveiter

  !----------------------------------------------------------------------------!
  ! write some global information

  call wpar_gl(upar)

  write (upar, fmpari) 'niter', nitert

  write (upar, fmhdr)  "disk information"
  write (upar, fmparfc) "alpha", alpha, "alpha parameter"
  if (cfg_magnetic) then
    write (upar, fmparfc) "eta", eta, "field rise parameter"
    write (upar, fmparfc) "eta", eta, "field rise parameter"
    write (upar, fmparfc) "nu", nu, "reconnection parameter"
    write (upar, fmparfc) "qcor", qcor, "corona parameter"
    write (upar, fmparfc) "qcor_fact", maxval(yy(c_qcor,:)), "same incl. pressure"
    write (upar, fmparfc) "xcor", qcor, "corona parameter (old name)"
  end if

  write (upar, fmpare) "radius", radius
  write (upar, fmpare) "radius_cm", radius * rschw
  write (upar, fmpare) "zscale", zscale
  write (upar, fmpare) "facc", facc
  write (upar, fmpare) "omega", omega

  write (upar, fmpare) "teff", teff
  write (upar, fmparf) "teff_keV", teff * keV_in_kelvin

  write (upar, fmpare) "zflux", (cgs_kapes * facc) / (cgs_c * omega**2)

  write (upar, fmpare) "rho_0", yy(c_rho,1)
  write (upar, fmpare) "temp_0", yy(c_temp,1)

  if (cfg_magnetic) then
    write (upar,fmpare) 'beta_0', yy(c_pgas,1) / yy(c_pmag,1)
    write (upar,fmpare) 'betakin_0', (yy(c_pgas,1) + yy(c_prad,1)) / yy(c_pmag,1)
    write (upar,fmpare) 'betakin_0_min', (yy(c_pgas,1) + yy(c_prad,1)) / yy(c_pmagmri,1)
    write (upar,fmpare) 'betakin_0_i', (2 * eta + alpha * (nu - 1)) / alpha
    write (upar,fmpare) 'alpha_0_fix1', yy(c_qmri,1) * alpha
    write (upar,fmpare) 'alpha_0_fix2', (2 * eta + alpha * nu) &
    & / (1 + (yy(c_pgas,1) + yy(c_prad,1)) / min(yy(c_pmagmri,1), yy(c_pmag,1)))
  end if

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
  write (upar, fmparl) "has_magnetic", cfg_magnetic
  write (upar, fmparl) "use_quench_mri", use_quench_mri
  write (upar, fmparl) "has_conduction", use_conduction
  write (upar, fmparl) "use_prad_in_alpha", use_prad_in_alpha
  write (upar, fmparl) "use_relcompt", use_relcompt
  write (upar, fmparl) "use_klein_nishina", use_klein_nishina

  !----------------------------------------------------------------------------!
  ! save the column information for col2python

  write_columns: block
    integer :: i

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
    real(dp) :: zphot, ztherm, zeqbc, ztmin, zinstabil

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

    ! ! vertical scale - weighted by gas pressure
    ! diskscale = sqrt(integrate(yy(c_pgas,:) * x**2, x) &
    !               /  integrate(yy(c_pgas,:),        x))
    ! write (upar, fmparec) 'zdisk_pgas', diskscale, &
    !      'disk height weighted by gas pressure'
    ! write (upar, fmparf) 'hdisk_pgas', diskscale / zscale
    
    if (cfg_magnetic) then
      ! vertical scale - weighted by magnetic pressure
      diskscale = sqrt(integrate(yy(c_pmag,:) * x**2, x) &
                    /  integrate(yy(c_pmag,:),        x))
      write (upar, fmparec) 'zdisk_pmag', diskscale, &
           'disk height weighted by magnetic pressure'
      write (upar, fmparf) 'hdisk_pmag', diskscale / zscale
    end if
    
    ! ! vertical scale - weighted by radiation pressure
    ! diskscale = sqrt(integrate(yy(c_prad,:) * x**2, x) &
    !               /  integrate(yy(c_prad,:),        x))
    ! write (upar, fmparec) 'zdisk_prad', diskscale, &
    !      'disk height weighted by radiation pressure'
    ! write (upar, fmparf) 'hdisk_prad', diskscale / zscale
    !
    ! write (upar, fmparf) 'hheat', integrate(yy(c_heat,:) * x, x) &
    !     / integrate(yy(c_heat,:), x) / zscale

    !--------------------------------------------------------------------------!

    call save_interpolated(0.0_dp, 'midpl', 'midplane')

    zphot = -1
    call tabzero(x, yy(c_tau,:), 1.0_dp, zphot)
    call save_interpolated(zphot, 'phot', 'tau = 1')

    ztherm = -1
    call tabzero(x, yy(c_tauth,:), 1.0_dp, ztherm)
    call save_interpolated(ztherm, 'therm', 'tau* = 1')

    if ( has_corona ) then
      zeqbc = -1
      call tabzero(x(ngrid:1:-1), yy(c_compfr,ngrid:1:-1), 0.5_dp, zeqbc)
      call save_interpolated(zeqbc, 'eqbc', 'compt. == brehms.')

      ztmin = -1
      call findtempmin(x, yy(c_temp,:), ztmin)
      call save_interpolated(ztmin, 'tmin', 'temperature minimum')

      call save_interpolated(max(ztherm, ztmin), 'cor', 'max(ztherm, ztmin)')
      call save_interpolated(max(ztherm, zeqbc), 'cor2', 'max(ztherm, zeqbc)')
      call save_interpolated(max(ztherm, zeqbc, ztmin), 'cor3', 'max(ztherm, zeqbc, ztmin)')
      
      zinstabil = -1
      if (sum(yy(c_instabil, ngrid-3:ngrid)) > 0) then
        call tabzero(x(ngrid:2:-1), yy(c_instabil, ngrid:2:-1), 0.0_dp, zinstabil)
        call save_interpolated(zinstabil, 'instabil', 'instability top')
        call save_interpolated(max(ztherm, zinstabil), 'therm_c', 'max(ztherm, zinstabil)')
        call save_interpolated(max(ztherm, zinstabil), 'tmin_c',  'max(ztmin,  zinstabil)')
        call save_interpolated(max(ztherm, ztmin, zinstabil), 'cor_c', 'max(ztherm, ztmin, zinstabil)')
        call save_interpolated(max(ztherm, zeqbc, zinstabil), 'cor2_c', 'max(ztherm, zeqbc, zinstabil)')
        call save_interpolated(max(ztherm, zeqbc, ztmin, zinstabil), 'cor3_c', 'max(ztherm, zeqbc, ztmin, zinstabil)')
      end if
    end if

    !--------------------------------------------------------------------------!

    write (upar, fmhdr)  "some global parameters"
    write(upar, fmparg) 'instabil', minval(yy(c_instabil,:))
    write(upar, fmparg) 'instabv', minval(yy(c_instabv,:))
    write (upar, fmparfc) 'compy_0', yy(c_compy,1), 'total Y'
    write (upar, fmpare) 'frad_top', yy(c_frad, ngrid)
    write (upar, fmpare) 'fmag_top', yy(c_fmag, ngrid)
    write (upar, fmparf) 'fbfrac_top', yy(c_fbfr, ngrid)

    !--------------------------------------------------------------------------!

    write_max_rhograd: block
      use slf_deriv, only: loggrad
      real(dp) :: rhograd(ngrid)
      call loggrad(x, yy(c_rho,:), rhograd)
      write (upar, fmparfc) 'rhograd_max', maxval(rhograd), 'maximum density gradient'
      write (upar, fmparfc) 'rhograd_min', minval(rhograd), 'minimum density gradient'
    end block write_max_rhograd

    if (cfg_magnetic) then
      write_quench: block
        real(dp) :: x1(ngrid), x2(ngrid), qmri_avg

        if (use_quench_mri) then
          x1(:) = yy(c_pgas, :) + yy(c_pmag, :) &
          &   + merge(yy(c_prad, :), 0._dp, use_prad_in_alpha)
          x2(:) = yy(c_qmri, :) * x1(:)
          qmri_avg = integrate(x2, x) / integrate(x1, x)
        else
          qmri_avg = 1
        end if

        write(upar, fmparf) 'qmri_avg', qmri_avg
        write(upar, fmparf) 'alphaeff', qmri_avg * alpha
      end block write_quench
    end if

  end block write_disk_globals

  !----------------------------------------------------------------------------!

  if (with_perf) then
    call system_clock(timing(2))
    print '("PERF", 1x, f7.4)', (timing(2) - timing(1)) / 1d9
  end if

  !----------------------------------------------------------------------------!
  ! clean up

  close(upar)
  if (.not. converged) error stop "not converged :("

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

    if (z < 0 .or. .not. ieee_is_normal(z)) return

    write (upar, fmhdr) 'properties in ' // comment

    write (upar, fmparec) 'z' // keyword, z, "location of" // comment
    write (upar, fmparf) 'h' // keyword, z / zscale

    call interpol(x, yy(c_temp,:), z, yz)
    write (upar,fmparec) 'temp_' // keyword, yz, 'temperature in ' // comment
    write (upar,fmparf) 'temp_' // keyword // '_keV', yz * keV_in_kelvin

    call interpol(x, yy(c_rho,:), z, yz)
    write (upar,fmparec) 'rho_' // keyword, yz, 'density in ' // comment
    write (upar,fmparec) 'rho_' // keyword // '_nh', yz / cgs_mhydr, &
      'density in ' // comment // ' (nH)'

    call interpol(x, yy(c_trad,:), z, yz)
    write (upar,fmparec) 'trad_' // keyword, yz, 'rad temp in ' // comment
    write (upar,fmparf) 'trad_' // keyword // '_keV', yz * keV_in_kelvin

    ! average temperature of the corona
    yz = interpolf(x, yy(c_tavg,:), z)
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
    write (upar,fmparec) 'beta_' // keyword, interpolf(x, yy(c_beta,:), z), &
        & 'magnetic beta in ' // comment
    write (upar,fmpare) 'ionxi_' // keyword, interpolf(x, yy(c_ionxi,:), z)

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

    call interpol(x, yy(c_nhtot,:), z, yz)
    write (upar, fmparec) 'nhtot_' // keyword, yz, 'column density above'

    write (upar, fmparg) 'kapsct_' // keyword, interpolf(x, yy(c_ksct,:), z)
    write (upar, fmparg) 'kapabs_' // keyword, interpolf(x, yy(c_kabs,:), z)
    write (upar, fmparg) 'kapabp_' // keyword, interpolf(x, yy(c_kabp,:), z)
    write (upar, fmpargc) 'mag_gauss_' // keyword, &
    &  sqrt(8 * pi * interpolf(x, yy(c_pmag,:), z)), 'field strength in Gauss'

    write (upar, fmparg) 'pgas_' // keyword, interpolf(x, yy(c_pgas,:), z)
    write (upar, fmparg) 'prad_' // keyword, interpolf(x, yy(c_prad,:), z)
    write (upar, fmparg) 'pmag_' // keyword, interpolf(x, yy(c_pmag,:), z)
    write (upar, fmparg) 'pmagmri_' // keyword, interpolf(x, yy(c_pmagmri,:), z)


  end subroutine

  !----------------------------------------------------------------------------!
  ! save each relaxation iteration to a numbered file (001, 002, 003 etc)

  subroutine saveiter(iter)

    integer, intent(in), optional :: iter
    character(len = 256) :: fn
    integer :: i

    if (present(iter)) then
      write(fn, '(A,".",I0.3,".dat")') trim(outfn), iter
    else 
      fn = trim(outfn) // '.dat'
    endif

    call fillcols(yv, c_, yy)
    
    open(33, file = trim(fn), action = 'write')
    write(33, '("#", a5, 2a14, *(a14))') 'N', 'Z', 'H', (trim(labels(i)), i = 1, size(labels))

    do i = 1, ngrid
      write(33, '(i6, *(es14.5e3))') i, x(i), x(i) / zscale, yy(:,i)
    end do

    close(33)

  end subroutine

  !----------------------------------------------------------------------------!

  subroutine evaluate_fout(model, x, yv, yy)
    integer, intent(in) :: model
    real(dp), intent(in) :: x(:), yv(:,:)
    real(dp), intent(inout) :: yy(:,:)
    procedure(funout_t), pointer :: fout
    integer :: i, n

    n = size(x)

    ! select the output function appropriate for this model
    call mrx_sel_fout(model, fout)

    ! evaluate the function for each point
    do concurrent (i = 1:n)
      call fout(x(i), yv(:,i), yy(:,i))
    end do
  end subroutine evaluate_fout

  !----------------------------------------------------------------------------!

  ! subroutine post_corona(model, x, yv, yy, t)
  !   use heatbalance, only: heatbil2
  !   integer, intent(in) :: model
  !   real(dp), intent(in) :: x(:)
  !   real(dp), intent(inout) :: yv(:,:), yy(:,:)
  !   real(dp) :: temp_old, t
  !   integer :: c_(6)
  !   integer :: i

  !   call mrx_sel_hash(model, c_)

  !   do concurrent (i = 1:ngrid)
  !     temp_old = yy(c_temp,i)
  !     call heatbil2(yy(c_rho,i), yy(c_temp,i), yy(c_trad,i), &
  !           yy(c_heat,i), .false.)
  !     if (t /= 1) yy(c_temp,i) = yy(c_temp,i)**t * temp_old**(1-t)
  !   end do

  !   yv(c_(2),:) = yy(c_temp,:)

  !   call evaluate_fout(model, x, yv, yy)

  ! end subroutine post_corona

  !----------------------------------------------------------------------------!

  subroutine fillcols(yv,c_,yy)
    use slf_deriv, only: loggrad
    use heatbalance, only: heatbil2
    real(dp) :: kabs,ksct,kabp,rhom,tempm,tradm,tcorrm,dx
    real(dp), dimension(:,:), intent(in) :: yv
    integer, dimension(:), intent(in) :: c_
    real(dp), dimension(:,:), intent(inout) :: yy
    integer :: i

    call evaluate_fout(model, x, yv, yy)

    if (cfg_post_corona .and. cfg_simple_post) then
              ! in simple mode, only do one iteration and correct gas pressure
        ! otherwise heating might get broken
      do concurrent (i = 1:ngrid)
        call heatbil2(yy(c_rho,i), yy(c_temp,i), yy(c_trad,i), &
              yy(c_heat,i), .false.)
      end do
      yy(c_pgas,:) = yy(c_pgas,:) * yy(c_temp,:) / yy(c_trad,:)

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)
    end if

    yy(c_dnh,:) = yy(c_rho,:) / cgs_mhydr

    ! opacities
    yy(c_ksct,:) = fkesp(yy(c_rho,:), yy(c_temp,:))
    yy(c_kabs,:) = fkabs(yy(c_rho,:), yy(c_temp,:))
    yy(c_kabp,:) = fkabp(yy(c_rho,:), yy(c_temp,:))
    yy(c_kcnd,:) = fkcnd(yy(c_rho,:), yy(c_temp,:))

    ! cooling components: brehmstrahlung and compton
    cooling: block
      yy(c_coolb,:) = 4 * cgs_stef * yy(c_rho,:) * yy(c_kabp,:)   &
            * yy(c_temp,:)**4

      yy(c_coolc,:) = 4 * cgs_stef * yy(c_rho,:) * yy(c_ksct,:)   &
          * yy(c_trad,:)**4 * cgs_k_over_mec2 * 4 * yy(c_temp,:) &
          * merge(1 + 4 * cgs_k_over_mec2 * yy(c_temp,:), &
          & 1.0_dp, use_relcompt)

      yy(c_compfr,:) = yy(c_coolc,:) / (yy(c_coolb,:) + yy(c_coolc,:))
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

    integrate_tau: do i = ngrid-1, 1, -1
      rhom = (yy(c_rho,i) + yy(c_rho,i+1)) / 2
      if (rhom < 0) rhom = 0
      tempm = (yy(c_temp,i) + yy(c_temp,i+1)) / 2
      tradm = (yy(c_trad,i) + yy(c_trad,i+1)) / 2
      dx = x(i+1) - x(i)

      kabs = merge(fkabs(rhom,tempm), 0.0_dp, tempm > 0)
      kabp = merge(fkabp(rhom,tempm), 0.0_dp, tempm > 0)
      ksct = merge(fkesp(rhom,tempm), 0.0_dp, tempm > 0)

      yy(c_tau,  i) = yy(c_tau,  i+1) + dx * rhom * kabs
      yy(c_taues,i) = yy(c_taues,i+1) + dx * rhom * ksct
      yy(c_tauth,i) = yy(c_tauth,i+1) + dx * rhom * sqrt(kabp * (kabp + ksct))

      yy(c_tavg, i) = yy(c_tavg, i+1) + dx * rhom * kabs * tempm

      tcorrm = merge(1 + 4 * cgs_k_over_mec2 * tempm, 1.0_dp, &
        use_relcompt)

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
    yy(c_tavg,ngrid) = yy(c_temp,ngrid)
    yy(c_tavg,:ngrid-1) = yy(c_tavg,:ngrid-1) / yy(c_tau,:ngrid-1)

    yy(c_radpz,:) = yy(c_kabs,:) * yy(c_frad,:) / (cgs_c * omega**2 * x)
    if (x(1) == 0) yy(c_radpz,1) = yy(c_radpz,2)

    ! magnetic beta parameter
    if (cfg_magnetic) then
      yy(c_beta,:) = yy(c_pgas,:) / yy(c_pmag,:)
      yy(c_betagen,:) = (yy(c_pgas,:) + merge(yy(c_prad,:), 0._dp, use_prad_in_alpha)) / yy(c_pmag,:)
      ! yy(c_betamri,:) = sqrt(yy(c_pgas,:) / (1.67 * yy(c_rho,:))) * 2 / (omega * radius * rschw)
      yy(c_pmagmri,:) = sqrt(1.67 * yy(c_pgas,:) * yy(c_rho,:)) * (omega * radius * rschw) / 2
    end if

    ! here we compute d ln (cooling) / d ln T
    instability: block
      use heatbalance, only: fcool2
      real(dp) :: ptot(ngrid)
      call fcool2(yy(c_rho,:), yy(c_temp,:), yy(c_trad,:), yy(c_cool,:), yy(c_cool_dr,:), yy(c_cool_dT,:))
      yy(c_cool_dr,:) = yy(c_rho,:)  * yy(c_cool_dr,:) / yy(c_cool,:)
      yy(c_cool_dT,:) = yy(c_temp,:) * yy(c_cool_dT,:) / yy(c_cool,:)
      yy(c_instabil,:) = yy(c_cool_dT,:) - yy(c_cool_dr,:)
      ptot(:) = yy(c_pgas,:) + yy(c_pmag,:) + yy(c_prad,:)
      yy(c_instabv,:) = yy(c_cool_dT,:) - (1 - yy(c_pmag,:) / ptot) * yy(c_cool_dr,:)
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

    if (cfg_magnetic) then
      call loggrad(x, yy(c_pmag,:), yy(c_qcor,:))
      yy(c_qcor,:) = -yy(c_qcor,:)
    end if

    ! column density
    yy(c_coldens,1) = 0
    integrate_coldens: do i = 1, ngrid-1
      yy(c_coldens, i+1) = yy(c_coldens,i) &
      + (yy(c_rho,i) + yy(c_rho,i+1)) * (x(i+1) - x(i)) / 2
    end do integrate_coldens

    yy(c_nhtot,ngrid) = 0
    integrate_nhtot: do i = ngrid-1, 1, -1
      yy(c_nhtot, i) = yy(c_nhtot,i+1) &
      + (yy(c_rho,i) + yy(c_rho,i+1)) * (x(i+1) - x(i)) / (2 * cgs_mhydr)
    end do integrate_nhtot

    instabil_ref: block
      real(dp) :: wsp1, wsp2

      wsp1 = (8. / 3.) * cgs_kapes / kram0p(abuX, abuZ) * cgs_k_over_mec2 * 1e27_r64
      wsp2 = 4 * cgs_stef * (cgs_kapes * 4 * cgs_k_over_mec2)**2 / kram0p(abuX, abuZ) * 1e57_r64

      yy(c_rho_max,:) = wsp1 * (yy(c_trad,:) / 1e6)**5 / (yy(c_temp,:) / 1e6)**0.5
      yy(c_cool_ref,:) = wsp2 * (yy(c_trad,:) / 1e6)**9.5
      yy(c_cf_x,:) = yy(c_rho,:) / (wsp1 * (yy(c_trad,:) / 1e6)**4.5)
      yy(c_cf_y,:) = yy(c_temp,:) / yy(c_trad,:)
    end block instabil_ref

  end subroutine

  !----------------------------------------------------------------------------!
  ! searches for temperature minimum in a given array

  pure subroutine findtempmin(x,temp,xtmin)
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

  subroutine generate_grid(x, ztop)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: ztop
    integer :: i, ngrid
    real(dp) :: z0

    ngrid = size(x)
    z0 = ztop * grid_nonlnr

    select case (tgrid)
    case (grid_linear)
      do i = 1, ngrid
        x(i)  = space_linear(i, ngrid, ztop)
      end do
    case (grid_log)
      do i = 1, ngrid
        x(i)  = space_linlog(i, ngrid, ztop / z0) * z0
      end do
    case (grid_linlog)
      call space_linlog2(x(:), ztop / z0)
      x(:) = x * z0
    case (grid_asinh)
      do i = 1, ngrid
        x(i)  = space_asinh(i, ngrid, ztop / z0) * z0
      end do
    case (grid_pow2)
      do i = 1, ngrid
        x(i)  = space_pow2(i, ngrid, ztop / z0) * z0
      end do
    case default
      error stop "this grid is not supported"
    end select

  end subroutine generate_grid

  !----------------------------------------------------------------------------!

  subroutine rdconf(cfg)
    integer :: errno
    type(config), intent(inout) :: cfg
    character(len=2048) :: buf

    call mincf_get(cfg, "alpha", buf, errno)
    if ( iand(errno, mincf_not_found) .eq. 0 )  then
      read (buf,*) alpha
    else
      error stop "alpha-parameter (key: alpha) is REQUIRED!"
    end if

    call mincf_get(cfg, "radius", buf, errno)
    if ( iand(errno, mincf_not_found) .eq. 0 )  then
      read (buf,*) radius
    else
      error stop "Radius (key: radius) is REQUIRED!"
    end if

    call mincf_get(cfg, "nu", buf, errno)
    if ( iand(errno, mincf_not_found) .eq. 0 )  then
      read (buf,*) nu
    else
      nu = 0.1 / alpha
      write (0, '("no nu given, assuming nu = ",f5.2)') nu
    end if

    call mincf_get(cfg, "eta", buf, errno)
    if ( iand(errno, mincf_not_found) .eq. 0 )  then
      read (buf,*) eta
      if (eta < 0) error stop 'eta < 0'
    else
      eta = 0.5 * (alpha / 0.5)**0.5
      write (0, '("no eta given, assuming alpha/eta = ",f5.3,"/",f5.3)') alpha, eta
    end if

  end subroutine

  !----------------------------------------------------------------------------!

  elemental function err_ramp(x, y1, p) result(y)
    real(dp), intent(in) :: x, y1, p  
    real(dp) :: a, y

    a = (1 / y1)**(1 / p) - 1
    y = 1 / (1 + a * x)**p
  end function

  !----------------------------------------------------------------------------!

  subroutine rdargvrx
    integer :: i
    character(2**8) :: arg

    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)

      ! write every iteration to another file? useful to demonstrate relaxation
      case ("-write-all", "-all")
        cfg_write_all_iters = .TRUE.

      ! recalculate the cooling-heating balance after relaxation? works best
      ! with -compton switch or alone (not much sense with -corona switch)
      case ("-post-corona", "-post")
        cfg_post_corona = .TRUE.
        cfg_simple_post = .false.
      case ("-no-post-corona", "-no-post")
        cfg_post_corona = .FALSE.
      case ("-simple-post-corona", "-simple-post")
        cfg_post_corona = .TRUE.
        cfg_simple_post = .true.
      case ("-no-simple-post-corona", "-no-simple-post")
        cfg_post_corona = .FALSE.

      case ("-estim-trad", "-no-estim")
        cfg_corestim_method = "none"
      case ("-estim-new", "-estim-full")
        cfg_corestim_method = "full"
      case ("-estim-instab")
        cfg_corestim_method = "instab"
      case ("-estim-old", "-estim-compt")
        cfg_corestim_method = "compt"

      ! enable PP condition for MRI shutdown? can cause trouble for convergence
      case ("-quench", "-qmri")
        use_quench_mri = .true.
      case ("-no-quench", "-no-qmri")
        use_quench_mri = .false.

      ! use P_rad in alpha prescription?
      case ("-prad-alpha", "-alpha-prad")
        use_prad_in_alpha = .TRUE.
      case ("-no-prad-alpha", "-no-alpha-prad")
        use_prad_in_alpha = .FALSE.

      ! use flux correction for highly magnetized disks?
      case ("-fluxcorr", "-flux-correction")
        use_flux_correction = .TRUE.
      case ("-no-fluxcorr", "-no-flux-correction")
        use_flux_correction = .FALSE.

      case('-alpha')
        cfg_magnetic = .false.

      case ("-perf","-with-perf")
        with_perf = .true.

      case ('-trim-vacuum', '-trim')
        cfg_trim_vacuum = .true.
      case ('-no-trim-vacuum', '-no-trim')
        cfg_trim_vacuum = .false.

      case ('-smooth')
        cfg_smooth_delta = .true.
      case ('-no-smooth')
        cfg_smooth_delta = .false.

      end select
    end do
  end subroutine

  !----------------------------------------------------------------------------!

  pure subroutine smooth2(x, nit)
    real(dp), intent(inout) :: x(:)
    real(dp) :: x1(size(x))
    integer :: n, i, it
    integer, intent(in) :: nit

    n = size(x)
    
    do it = 1, nit
      x1(1:n) = x(1:n)
      do concurrent(i = 2:n-1)
        x(i) = (x1(i-1) + 2 * x1(i) + x1(i+1)) / 4
      end do
    end do

  end subroutine

  !----------------------------------------------------------------------------!

end program
