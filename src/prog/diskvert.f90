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
  logical :: user_ff, user_bf, converged = .false., has_corona
  integer, dimension(6) :: c_
  integer, parameter :: upar = 92
  !----------------------------------------------------------------------------!
  integer(int64) :: timing(2)
  logical :: with_perf = .false.
  !----------------------------------------------------------------------------!
  logical :: cfg_write_all_iters = .FALSE., cfg_write_corona_only = .false.
  character, parameter :: EQUATION_SIMPBALANCE = 'D'
  logical :: cfg_smooth_delta = .false.
  logical :: cfg_post_corona = .false., cfg_iter_post = .false., &
    cfg_magnetic = .true., cfg_trim_vacuum = .false.
  character(len=8) :: cfg_corestim_method = ""
  !----------------------------------------------------------------------------!
  real(dp), parameter :: trim_density_thresh = 1e-11
  real(dp), parameter :: typical_hdisk = 12
  real(dp) :: logstep = -1, hydrox_max = 1
  !----------------------------------------------------------------------------!


  integer, parameter :: ncols  = n_yout + 39, &
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
      c_instabp   = n_yout + 21, &
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
      c_instabm   = n_yout + 35, &
      c_coolerr   = n_yout + 36, &
      c_dr_pgas   = n_yout + 37, &
      c_dr_prad   = n_yout + 38, &
      c_dr_pmag   = n_yout + 39

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
  labels(c_instabp) = 'instabp'
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
  labels(c_instabm) = 'instabm'
  labels(c_coolerr) = 'coolerr'
  labels(c_dr_pgas) = 'dr_pgas'
  labels(c_dr_prad) = 'dr_prad'
  labels(c_dr_pmag) = 'dr_pmag'

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

    if (qcor <= 1 .and. .not. use_quench_mri) write (*, *) 'warning: qcor <= 1'
    ! if (qcor <= 1 .and. use_flux_correction .and. .not. use_quench_mri) &
    ! &   error stop "qcor <= 1 .and. use_flux_correction .and. .not. use_quench_mri"
    if (qcor <= 0) error stop "qcor <= 0"
  end if

  !----------------------------------------------------------------------------!
  ! estimate the interval height

  if (cfg_auto_htop) then
    block_auto_height: block
      real(dp) :: hg, hr, hm, hm1, qcor1
      real(dp), parameter :: cutoff = 1e-7

      
      hg = zdisk_ss73 / zscale
      hr = cgs_kapes * facc / (omega**2 * cgs_c * zscale)
      write(*, '("SS73 height    ", f10.1)') hg
      write(*, '("rad press heigh", f10.1)') hr
      
      if (cfg_magnetic) then
        qcor1 = (2 * eta + alpha * (nu - merge(1e-2_dp, 1.0_dp, use_quench_mri))) / eta
        hm1 = sqrt((4 + alpha * nu / eta) * (cutoff**(-2 / (qcor1 + 1)) - 1))
        hm = 1 + 1 / cutoff**(1 / (qcor1 + 2.5_dp))
        write(*, '("magnetic height", f10.1)') hm, hm1
        htop = 12 * hg * hm
      else
        htop = 9 * hg
      end if

      htop = htop * merge(1.1, 1.0, cfg_trim_vacuum)
      
      write(*, '("assumed height ", f10.1)') htop
    end block block_auto_height
  end if
  write(upar, fmpare) 'htop_init', htop
    
  !----------------------------------------------------------------------------!

  if (logstep < 0) then
    if (cfg_magnetic) then
      logstep = logstep_q(qcor, -3.1_dp, 0.7_dp)
    else
      logstep = 1e-3
    end if
  end if

  if (ngrid < 300) error stop 'ngrid < 300'
  if (ngrid > 4000) error stop 'ngrid > 4000'
  ngrid = 16 * nint(ngrid / 16.)

  !----------------------------------------------------------------------------!

  ! get the model number
  model = mrx_number( 'D', cfg_magnetic, use_conduction )
  call mrx_sel_nvar(model, ny)
  call mrx_sel_hash(model, C_)

  allocate(x(ngrid), Y(ny*ngrid), YY(ncols,ngrid))

  !----------------------------------------------------------------------------!
  ! generate the grid

  call generate_grid(tgrid, x, htop * zscale, logstep)

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
        y_pmag(i) = max(y_pmag(i), 2 * pth / (beta_0 + nu))
      end if
    end do

    if (cfg_write_all_iters .and. (.not. cfg_write_corona_only)) call saveiter(0)

  end block initial_profile

  !----------------------------------------------------------------------------!

  relaxation_block: block

    real(dp), allocatable, target :: dY(:), M(:,:)
    real(dp), allocatable :: MB(:,:)
    integer, dimension(:), allocatable :: ipiv
    logical, dimension(:), allocatable :: errmask
    real(dp) :: err, err0, ramp, rholim, ers, rho_floor
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
    cndredu = 1

    converged = .false.

    write(*, '(a5, 1x, a10, 10x, 3(2x, a5, "%"), 2x, a5)') &
    & 'ITR', 'ERR', 'RAMP', 'R_OPA', 'R_MRI', 'NGRHO'

    relx_disk : do iter = 1, 2222

      ! negative rho is number of recent iterations where rho was negative or posivie
      if (any(y_rho <= 0)) then
        negative_rho = max(0, negative_rho) + 1
      else
        negative_rho = merge(negative_rho - 1, 0, negative_rho <= 0 .and. iter > 1)
      end if

      if ((opacities_kill <= 0 .and. iter > 5 .and. err < 0.03 .and. negative_rho < -3) &
      &   .or. (opacities_kill > 0 .and. err < 0.3)) then
        opacities_kill = min(1._dp, (opacities_kill + 0.3) * 1.05 - 0.3)
      end if

      if (iter > 500 .and. err < 1e-5 .and. negative_rho > 0) exit relx_disk

      call mrx_matrix(model, x, Y, M, dY)
      call mrx_bandim(model, kl, ku)
      call m2band(M, KL, KU, MB)
      call dgbsv(size(MB,2), KL, KU, 1, MB, size(MB,1), ipiv, dY, size(dY), errno)

      errmask(:) = (Y .ne. 0) .and. ieee_is_normal(dY)
      err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))

      if (ieee_is_nan(err) .or. (iter > 10 .and. (err > 1e10))) then
        write(*, '("'// achar(27) //'[1;31;7mdiverged: ", Es9.2, " -> ", Es9.2, "'&
            // achar(27) //'[0m")') err0, err
        converged = .false.
        exit relaxation_block
      end if

      ers = merge((err / err0)**0.33, 1._dp, err0 > 0)
      ers = max(0.33_dp, min(8.0_dp, ers)) * merge(1._dp, 3._dp, negative_rho < -2)
      ramp = max(err_ramp(err * ers, 0.25_dp, 0.50_dp), 1e-3_dp)

      ! rholim = minval(-yv(:,c_(1)) / dyv(:,c_(1)), &
      ! & mask=yv(:,c_(1)) > 0 .and. dyv(:,c_(1)) < 0 .and. ieee_is_normal(-yv(:,c_(1)) / dyv(:,c_(1))))
      ! rholim = min(1._dp, rholim)
      ! if (rholim < 1) ramp = ramp * min(1._dp, max(0.1_dp, 0.8 * rholim / ramp))
        
      write(*, '(a, i5, 1x, 1es10.2, 4x, "ramps=", 3(2x, F5.1, "%"), 2x, i5, 1x, a, a)') &
      & trim(merge(achar(27) // '[33;1m', repeat(' ', 7), err0 < err .and. err0 > 0)), &
      & nitert+1, err,  100*ramp, 100*opacities_kill, 100*qmri_kill, negative_rho,  &
      & trim(merge('rho<0', repeat(' ', 5), negative_rho > 0)), &
      & trim(merge(achar(27) // '[0m', repeat(' ', 4), err0 < err .and. err0 > 0))

      Y(:) = Y + dY * ramp

      ! fix negative densities
      if (negative_rho > 5) then
        where (y_rho < 0) y_rho = -y_rho
        call smooth2(y_rho, 8)
      end if

      nitert = nitert + 1
      if (cfg_write_all_iters .and. (.not. cfg_write_corona_only)) call saveiter(nitert)

      if (err < 0.03 .and. use_quench_mri .and. negative_rho < -2) then
        qmri_kill = min(1._dp, (qmri_kill + 0.5) * 1.03 - 0.5)
        threshpow =  1 + 3 * qmri_kill
      end if

      if (err < 1e-6 .and. (err0 / err > 5 .or. err < 1e-9) .and. opacities_kill > 0.999 &
      &   .and. (qmri_kill > 0.999 .or. .not. use_quench_mri) .and. negative_rho < 0) then
        write(*, '("STRUCTURE convergence reached with error = '// achar(27) &
            // '[1m",ES9.2,"'// achar(27) //'[0m")') err
        converged = .true.
        exit relx_disk
      end if

      err0 = err

      rho_floor = trim_density_thresh * maxval(y_rho)

      if (cfg_trim_vacuum .and. (iter > 5 .and. opacities_kill > 0.3) &
      & .and. (negative_rho < -9) &
      & .and. any(y_rho(ngrid/2:) < rho_floor)) then
        trim_space_vacuum: block
          real(dp) :: xcut, xnew(ngrid)

          call tabzero(x, y_rho, 1.5 * rho_floor, xcut)
          if (xcut / zscale < 3 .or. xcut / x(ngrid) > 0.99) exit trim_space_vacuum

          write(*, '("'// achar(27) //'[96;1m<<< trimming top from ", &
          &   f0.1, " to ", f0.1, "'// achar(27) //'[0m")') &
          &   x(ngrid) / zscale, xcut / zscale

          call generate_grid(tgrid, xnew, xcut, logstep)
          call interpolate_grid(xnew, x, y, ny)
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
      
      
      set_new_grid: block
        real(dp), allocatable :: xnew(:)
        integer :: ngrid_new

        ngrid_new = 16 * nint((ngrid * 1.33) / 16)

        allocate(xnew(ngrid_new))

        ! 0.1: 28/111, 0.5: 26/111, 0.01: 30, auto: 24
        call generate_grid(grid_log, xnew, x(ngrid), logstep_q(qcor, -1.7_dp, 2.4_dp))
        call interpolate_grid(xnew, x, y, ny)

        ngrid = ngrid_new
      end block set_new_grid
 

      deallocate(dY, M, errmask, ipiv, yy)
      allocate(dY(ny*ngrid), M(ny*ngrid,ny*ngrid), errmask(ny*ngrid), &
          ipiv(ny*ngrid), yy(ncols, ngrid))
      M(:,:) = 0
      
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
        real(dp) :: heat(ngrid), tmax, wsp1, tmax2, tcompt
        real(dp) :: mixin = 0.5
        integer :: i


        call evaluate_fout(model, x, yv, yy)
        heat(:) = yy(c_heat, :)

        ! y_pgas(:) = cgs_k_over_mh / miu * y_rho(:) * y_trad(:)
        ! y_prad(:) = cgs_a * y_trad(:)**4 / 3
        ! heat(:) = 2 * (eta + alpha * nu) * omega * y_pmag(:) &
        !   - alpha * omega * (y_pgas(:) + y_pmag(:) &
        !   + merge(y_prad(:), 0.0_dp, use_prad_in_alpha))
        ! heat(:) = (2 * eta + alpha * (2 * nu - 1)) * omega * y_pmag(:)

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
          mixin = 0.6
          do concurrent (i = 1:ngrid)
            tcompt = heat(i) * (cgs_mel * cgs_c**2) &
            &   / (16 * cgs_boltz * (cgs_kapes * y_rho(i)) &
            &   * (cgs_stef * y_trad(i)**4) )
            y_temp(i) = y_trad(i) + tcompt
            ! y_temp(i) = sqrt(y_trad(i)**2 + (y_trad(i) + tcompt)**2)
          end do

        case ('full')
          mixin = 0.9
          do concurrent (i = 1:ngrid)
            call heatbil2(y_rho(i), y_temp(i), y_trad(i), heat(i), .false.)
          end do

        case ('instab')
          mixin = 0.9
          wsp1 = (8. / 3.) * cgs_kapes / kram0p(abuX, abuZ) * cgs_k_over_mec2
          do concurrent (i = 1:ngrid)
            call heatbil2(y_rho(i), y_temp(i), y_trad(i), heat(i), .false.)
            tmax = sqrt((1.38 * y_trad(i))**2 + ((wsp1 * y_trad(i)**5 / y_rho(i))**2)**2)
            ! tmax2 = 1.2 * (2 * eta / alpha + 2 * nu - 1) * y_pmag(i) / (cgs_k_over_mh * y_rho(i) / 0.5)
            y_temp(i) = min(y_temp(i), tmax)
            ! y_temp(i) = 1.06 / (1 / y_temp(i)**4 + 1 / tmax**4)**(1. / 4.) 
          end do
          
        case default
          y_temp(:) = y_trad(:)
          write(*, *) 'warning: no corona estimation'
        end select
        
        y_temp(:) = y_trad(:)**(1 - mixin) * y_temp(:)**mixin
      end block estimate_corona

      !--------------------------------------------------------------------------!

      nitert = nitert + 1
      if (cfg_write_all_iters) call saveiter(nitert)
      
      !--------------------------------------------------------------------------!
      
      converged = .false.
      err0 = 0
      ! cndredu = 0
      condux = 0.0

      hydrox = 1e-5

      relx_corona : do iter = 1, 999

        ! if (err > 1e-2) then
          ! block
          !   real(dp) :: a(ngrid)
          !   call random_number(a)
          !   ! y_rho(:) = y_rho * (1._dp + 1e-3_dp * (2 * a - 1))
          !   y_temp(:) = y_trad + (y_temp - y_trad) * (1._dp + 3e-6_dp * (2 * a - 1))
          ! end block
        ! end if

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
            call smooth2(dYv(i,:), 1)
          end do
        end if
        
        errmask(:) = (Y .ne. 0) .and. ieee_is_normal(dY)
        err = sqrt(sum((dY/Y)**2, errmask) / count(errmask))
        ers = merge((err / err0)**0.50, 1._dp, err0 > 0)
        ers = max(0.50_dp, min(8.0_dp, ers)) 
        ramp = max(err_ramp(err * ers, 0.10_dp, 0.50_dp), 1e-4_dp)
        
        ! (iter > 1 .and. err > err0)
        write(*, '(a, I5, 2X, ES9.2, 2X, F5.1, "%", 3x, f5.1, "%", a)') '', nitert+1, err, 100*ramp, 100*hydrox, ''

        if (iter > 1 .and. err < 0.5) then
          hydrox = min(hydrox_max, 1.2 * (hydrox + 0.25) - 0.25)
          condux = min(1.0_dp, 1.03 * (condux + 1e-2) - 1e-2)
          if (use_conduction)  print *, 'condux', condux
        end if
        
        where (.not. ieee_is_normal(dyv)) dyv = 0
        Y(:) = Y + dY * ramp
        ! where (.not. ieee_is_normal(Y)) Y = 0

        ! where (y_temp <= 0) y_temp = 1.1 * y_trad
        ! where (y_rho <= 0) y_rho = -y_rho
        ! !   if (any(y_rho <= 0)) then
        ! !   print *, 'oops'
        ! !   where (y_rho <= 0) y_rho = -y_rho
        ! !   call smooth2(y_rho, 1)
        ! ! end if
        
        ! if (err0 > 0 .and. err / err0 > 1.2) then
        !   call smooth2(y_temp, 1)
        !   call smooth2(y_rho, 1)
        ! end if

        nitert = nitert + 1
        if (cfg_write_all_iters) call saveiter(nitert)

        if (ieee_is_nan(err) .or. (iter > 1 .and. err > err0 * 1e8) .or. (err > 1e12)) then
          write(*, '("' // achar(27) // '[1;31;7mCORONA diverged: ", Es9.2, " -> ", Es9.2, "'&
              // achar(27) //'[0m")') err0, err
          if (iter > 2) write(*, '(a, i4)') 'ncorit = ', iter
          exit relaxation_block
        end if

        if (err < 5e-5 .and. (err0 / err > 8 .or. err < 1e-8) .and. (hydrox >= 0.999 * hydrox_max) &
        & .and. .not. (use_conduction .and. condux < 0.999)) then
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
    if (converged .and. cfg_post_corona .and. cfg_iter_post) then
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
          yv(c_(2),i) = yy(c_temp,i)**0.2 * temp_old**0.8
        end do

        err = sum((2 * (yv(c_(2),:) - yy(c_temp,:)) / (yv(c_(2),:) + yy(c_temp,:)))**2)
        write(*, '(a, 1x, i3, 1x, es8.1)') 'POST', iter, err
    
        call evaluate_fout(model, x, yv, yy)

        do concurrent (i = 1:ngrid)
          yy(c_heat,i) = max(0._dp, yy(c_heat,i))
        end do

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
    write (upar, fmparf) "etaxi", eta + alpha * nu / 2
    write (upar, fmparfc) "nu", nu, "reconnection parameter"
    write (upar, fmparfc) "qcor", qcor, "corona parameter"
    write (upar, fmparfc) "qcor_fact", maxval(yy(c_qcor,:)), "same incl. pressure"
  end if

  write (upar, fmpare) "radius", radius
  write (upar, fmpare) "radius_cm", radius * rschw
  write (upar, fmpare) "zscale", zscale
  write (upar, fmpare) "facc", facc
  write (upar, fmpare) "omega", omega
  write (upar, fmpare) 'htop', htop
  write (upar, fmparec) 'mdot_edd', mdot_edd(mbh), 'g/s'
  write (upar, fmparfc) 'mdot_cgs', log10(mdot * mdot_edd(mbh)), 'log g/s'

  write (upar, fmpare) "teff", teff
  write (upar, fmparf) "teff_keV", teff * keV_in_kelvin

  write (upar, fmpare) "zflux", (cgs_kapes * facc) / (cgs_c * omega**2)

  write (upar, fmpare) "rho_0", yy(c_rho,1)
  write (upar, fmpare) "temp_0", yy(c_temp,1)

  if (cfg_magnetic) then
    write (upar,fmpare) 'beta_0', yy(c_pgas,1) / yy(c_pmag,1)
    write (upar,fmpare) 'beta_0_fact', yy(c_pgas,1) / yy(c_pmag,1)
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

  ! 1 - kH / c Om2 > 0
  block
    use slf_deriv, only: loggrad
    real(dp) :: h0, densgrad(ngrid)
    h0 = yy(c_kabs,1) * yy(c_heat,1) / (cgs_c * omega**2)
    ! h1 = 3 * cgs_boltz * yy(c_rho,1) * cgs_c / (miu * cgs_mhydr * 16 * cgs_stef * yy(c_trad,1)**3)
    write (upar, fmparf) 'hole_0', 1 - h0
    
    call loggrad(x, yy(c_rho,:), densgrad)
    write (upar, fmpare) 'max_densgrad', maxval(densgrad, 1)
    write (upar, fmpare) 'min_densgrad', minval(densgrad, 1)
    ! write (upar, fmparf) 'hole_1', 1 - h0 * (1 + h1)
    ! write (upar, fmparf) 'hole_01', h1
  end block

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

  write(upar, fmparf) "rho_spread", log10(maxval(yy(c_rho,:)) / minval(yy(c_rho,:)))

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
    logical :: has_instability

    has_instability = any(yy(c_instabp,:) < 0)

    write(upar, fmparl) 'has_instability', has_instability

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

    block
      real(dp) :: zphot, ztherm, ztmin, zinstabil, zcor

      call save_interpolated(0.0_dp, 'midpl', 'midplane')

      zphot = -1
      call tabzero(x, yy(c_tau,:), 1.0_dp, zphot)
      call save_interpolated(zphot, 'phot', 'tau = 1')

      ztherm = -1
      call tabzero(x, yy(c_tauth,:), 1.0_dp, ztherm)
      call save_interpolated(ztherm, 'therm', 'tau* = 1')

      if ( has_corona ) then

        ztmin = -1
        call findtempmin(x, yy(c_temp,:), ztmin)
        call save_interpolated(ztmin, 'tmin', 'temperature minimum')
        
        zcor = max(ztherm, ztmin)
        call save_interpolated(zcor, 'cor', 'max(ztherm, ztmin)')

        ! call tabzero(x, yy(c_temp,:) - 1.3838 * yy(c_trad,:), 0._dp, zcor)
        ! zcor = max(zcor, ztherm)
        ! call save_interpolated(zcor, 'cor', 'T > 1.38*Trad')
        
        zinstabil = -1
        if (sum(yy(c_instabp, ngrid-1:ngrid)) > 0) then
          call tabzero(x(ngrid:2:-1), yy(c_instabp, ngrid:2:-1), 0.0_dp, zinstabil)
        else
          zinstabil = x(ngrid)
        end if
        
        call save_interpolated(zinstabil, 'instabp', 'instability top (Pgas)')
        call save_interpolated(max(zcor, zinstabil), 'cor_c', 'max(zcor, zinstabil)')
        call save_interpolated(merge((zcor + zinstabil) / 2, zcor, zinstabil > zcor), &
        &   'cor_a', '(zcor+zinstabil)/2')

        zinstabil = x(minloc(yy(c_instabp, :), 1))
        call save_interpolated(zinstabil, 'instabp_m', 'instability min (Pgas)')

        zinstabil = -1
        if (sum(yy(c_instabm, ngrid-3:ngrid)) > 0) then
          call tabzero(x(ngrid:2:-1), yy(c_instabm, ngrid:2:-1), 0.0_dp, zinstabil)
          call save_interpolated(zinstabil, 'instabm', 'instability top (Pgas+Pmag)')
        end if

      end if

    end block

    !--------------------------------------------------------------------------!

    if (has_instability) then
      instab_params: block
        logical :: instabil_mask(ngrid)
        real(r64), allocatable :: rho_inst(:), temp_inst(:), x_inst(:)

        instabil_mask(:) = yy(c_instabp,:) < 0
        if (count(instabil_mask) < 2) exit instab_params

        rho_inst = pack(yy(c_rho,:), instabil_mask)
        temp_inst = pack(yy(c_temp,:), instabil_mask)
        x_inst = pack(x(:), instabil_mask)

        write(upar, fmpare) 'instabp_tavg', integrate(temp_inst * rho_inst, x_inst) / integrate(rho_inst, x_inst)
        write(upar, fmparf) 'instabp_taues', maxval(yy(c_taues,:), instabil_mask) - minval(yy(c_taues,:), instabil_mask)
      end block instab_params
    end if

    !--------------------------------------------------------------------------!

    write (upar, fmhdr)  "some global parameters"
    write(upar, fmparg) 'instabp', minval(yy(c_instabp,:))
    write(upar, fmparg) 'instabm', minval(yy(c_instabm,:))
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
    character(len=*), intent(in) :: keyword, comment
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

    write (upar, fmparg) 'instabp_' // keyword, interpolf(x, yy(c_instabp,:), z)
    write (upar, fmparg) 'instabm_' // keyword, interpolf(x, yy(c_instabm,:), z)
    write (upar, fmparg) 'cool_dr_' // keyword, interpolf(x, yy(c_cool_dr,:), z)
    write (upar, fmparg) 'cool_dT_' // keyword, interpolf(x, yy(c_cool_dT,:), z)

    write (upar, fmparg) 'cf_x_' // keyword, interpolf(x, yy(c_cf_x,:), z)
    write (upar, fmparg) 'cf_y_' // keyword, interpolf(x, yy(c_cf_y,:), z)

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

    if (cfg_post_corona .and. (.not. cfg_iter_post)) then
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

      call fcool2(yy(c_rho,:), yy(c_temp,:), yy(c_trad,:), yy(c_cool,:), yy(c_cool_dr,:), yy(c_cool_dT,:))
      yy(c_cool_dr,:) = yy(c_rho,:)  * yy(c_cool_dr,:) / yy(c_cool,:)
      yy(c_cool_dT,:) = yy(c_temp,:) * yy(c_cool_dT,:) / yy(c_cool,:)

      yy(c_instabp,:) = yy(c_cool_dT,:) - yy(c_cool_dr,:)
      yy(c_instabm,:) = yy(c_cool_dT,:) - yy(c_cool_dr,:) / (1 + 2 * yy(c_pmag,:) / yy(c_pgas,:))

      yy(c_coolerr,:) = (yy(c_heat,:) - yy(c_cool,:)) / yy(c_heat,:)
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

    density_gradients: block
      use slf_deriv, only: diffx, fdiffx
      real(dp) :: dx2(ngrid)
    ! - 2 * Dx(d['pmag']) / (p.omega**2 * Dx(d['z']**2) * dv.cgs_mhydr)
      dx2(:) = - omega**2 * fdiffx(x**2) / 2
      yy(c_dr_pgas,:) = fdiffx(yy(c_pgas,:)) / dx2
      yy(c_dr_prad,:) = fdiffx(yy(c_prad,:)) / dx2
      yy(c_dr_pmag,:) = fdiffx(yy(c_pmag,:)) / dx2
    end block density_gradients

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

  subroutine generate_grid(tgrid, x, ztop, logstep)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: ztop, logstep
    integer, intent(in) :: tgrid
    integer :: i, ngrid
    real(dp) :: z0

    ngrid = size(x)
    z0 = ztop * logstep

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
      call space_linlog2(x(:), 1 + ztop / z0)
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

  subroutine interpolate_grid(xnew, x, y, ny)
    use slf_interpol, only: interpol

    real(dp) :: xcut
    real(dp), intent(in) :: xnew(:)
    real(dp), intent(inout), allocatable :: x(:), y(:)
    real(dp), allocatable ::  ycopy(:)
    integer, intent(in) :: ny
    integer :: i,j

    if (size(x) * ny /= size(y)) error stop

    ycopy = y(:)

    if (size(xnew) /= size(x)) then
      deallocate(y)
      allocate(y(size(xnew) * ny))
    end if

    do j = 1, ny
      associate (y_j => y(j::ny), ycopy_j => ycopy(j::ny))
        do i = 1, size(xnew)
          call interpol(x, ycopy_j, xnew(i), y_j(i))
        end do
      end associate
    end do

    x = xnew(:)
  end subroutine

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
    integer :: i, skip
    character(2**8) :: arg, buf

    skip = 0

    do i = 1, command_argument_count()

      if (skip > 0) then
        skip = skip - 1
        cycle
      end if

      call get_command_argument(i, arg)
      select case (arg)

      ! write every iteration to another file? useful to demonstrate relaxation
      case ("-all")
        cfg_write_all_iters = .true.
      case ("-corona-only")
        cfg_write_corona_only = .true.
      case ("-no-corona-only")
        cfg_write_corona_only = .false.

      ! recalculate the cooling-heating balance after relaxation? works best
      ! with -compton switch or alone (not much sense with -corona switch)
      case ("-iter-post-corona", "-iter-post")
        cfg_post_corona = .TRUE.
        cfg_iter_post = .true.
      case ("-no-iter-post-corona", "-no-iter-post")
        cfg_post_corona = .FALSE.
        cfg_iter_post = .false.
      case ("-post-corona", "-post")
        cfg_post_corona = .TRUE.
      case ("-no-post-corona", "-no-post")
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

      case ('-simple-hydro')
        hydrox_max = 1e-5
      case ('-no-simple-hydro')
        hydrox_max = 1.0

      case ('-logstep')
        call get_command_argument(i + 1, buf)
        read(buf, *) logstep
        skip = skip + 1

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

  ELEMENTAL FUNCTION LOGSTEP_Q(QCOR, A, S) RESULT(LOGSTEP)
    REAL(DP), INTENT(IN) :: QCOR, A, S
    REAL(DP) :: LOGSTEP
    LOGSTEP = log(qcor + 1._dp) / log(15._dp)
    LOGSTEP = max(0._dp, min(1._dp, LOGSTEP))
    LOGSTEP = 10**(A + S * (LOGSTEP - 0.5_dp))
  END FUNCTION

end program
