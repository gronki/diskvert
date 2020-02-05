program diskvert_ss73

  use globals
  use settings
  use fileunits
  use confort
  use iso_fortran_env, only: sp => real32, dp => real64, int64, stdout => output_unit
  use ss73solution, only: ss73_estim, ss73_refin

  implicit none
  type(config) :: cfg
  real(dp) :: alpha, facc, teff, zscale, zdisk_ss73, rho_0_ss73, temp_0_ss73, omega, radius, coldens
  integer :: upar
  
  outfn = ""
  call rdargvgl
  call mincf_read(cfg)
  call rdconfgl(cfg)
  call rdconf(cfg)
  call mincf_free(cfg)

  upar = stdout
  if (outfn /= "") open(newunit = upar, file = trim(outfn) // '.txt', action = 'write')

  kappa_abs_0 = kram0(abuX, abuZ)

  call cylinder(mbh, mdot, radius, rschw, omega, facc, teff, zscale)

  write(upar, fmparec) 'mbh', mbh, '* Msun'
  write(upar, fmparec) 'mbh_g', mbh * cgs_msun, '[g]'
  write(upar, fmpare) 'mdot', mdot
  write(upar, fmpare) 'radius', radius
  write(upar, fmpare) 'radius_cm', radius * rschw
  write(upar, fmpare) 'alpha', alpha
  write(upar, fmparec) 'zscale', zscale, '[cm]'
  write(upar, fmparec) 'teff', teff, '[K]'
  write(upar, fmpare) 'd', zscale / (radius * rschw)
  write(upar, fmpare) 'omega', omega
  write(upar, fmparec) 'vfi', omega * radius * rschw / 1e5, '[km / s]'
  write(upar, fmparfc) 'vfi_c', omega * radius * rschw / cgs_c, '* c'
  
  call ss73_estim(mbh, mdot, radius, alpha, rho_0_ss73, temp_0_ss73, zdisk_ss73)
  write(upar, fmparec) 'zdisk_init', zdisk_ss73, 'initial estimate'
  write(upar, fmparfc) 'hdisk_init', zdisk_ss73 / zscale, 'initial estimate'
  
  call ss73_refin(mbh, mdot, radius, alpha, rho_0_ss73, temp_0_ss73, zdisk_ss73)
  write(upar, fmparec) 'zdisk', zdisk_ss73, '[cm]'
  write(upar, fmparfc) 'hdisk', zdisk_ss73 / zscale, '* zscale'
  
  write(upar, fmparec) 'rho_midpl', rho_0_ss73, '[g / cm3] - midplane'
  write(upar, fmparec) 'temp_midpl', temp_0_ss73, '[K] - midplane'
  
  coldens = rho_0_ss73 * zdisk_ss73 * sqrt(pi) / 2
  write(upar, fmparec) 'coldens', coldens, '[g / cm2]'
  write(upar, fmparec) 'nhtot', coldens / cgs_mhydr, '[1 / cm2]'
  
  write(upar, fmparec) 'vufo', omega * zdisk_ss73 / 1e5, '[km / s] - estimate UFO speed'
  
  write(upar, fmpare) 'kappa_abs_0', kappa_abs_0
  write(upar, fmparl) 'converged', .true.
  if (upar /= stdout) close(unit = upar)
  
  contains
  
  subroutine rdconf(cfg)
    type(config), intent(inout) :: cfg
    character(len = 1024) :: buf
    integer :: errno

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

  end subroutine

end program