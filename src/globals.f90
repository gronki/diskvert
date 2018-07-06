module globals

  use iso_fortran_env, only: r64 => real64
  use slf_cgs

  implicit none !--------------------------------------------------------!

  real(r64), parameter, private :: X0 = 0.68d0 / cgs_kapes_hydrogen - 1
  ! real(r64), parameter, private :: X0 = 0.7381
  real(r64), parameter, private :: Z0 = 0.0134

  real(r64) :: abuX = X0
  real(r64) :: abuZ = Z0

  real(r64) :: kappa_es = cgs_kapes_hydrogen * (1 + X0) / 2
  real(r64) :: kappa_ff_0 = 3.68d22 * (1 - Z0) * (1 + X0)
  real(r64) :: kappa_bf_0 = 4.34d25 * Z0 * (1 + X0)
  real(r64) :: kappa_abs_0 = 3.68d22 * (1 - Z0) * (1 + X0)

  logical :: use_opacity_ff = .true.
  logical :: use_opacity_bf = .false.
  logical :: use_opacity_planck = .true.
  logical :: use_conduction = .false.

  real(r64) :: mbh, mdot, rschw

  real(r64) :: cndredu = 1e-1

  real(r64), parameter :: miu = 0.50
  real(r64), parameter :: pi = 4*atan(real(1,r64))

  !   którego równania użyc?
  character, parameter :: EQUATION_DIFFUSION = 'D', &
  EQUATION_BALANCE = 'C', EQUATION_COMPTON = 'W'
  character :: cfg_temperature_method = EQUATION_DIFFUSION

contains !-----------------------------------------------------------------!

  elemental subroutine cylinder(mbh, mdot, r, rschw, omega, facc, teff, zscale)
    real(r64), intent(in) :: mbh, mdot, r
    real(r64), intent(out) :: rschw, omega, facc, teff, zscale
    real(r64) :: GM, mdot_crit, mdot_edd

    GM = cgs_graw * mbh * sol_mass

    rschw = 2 * GM / cgs_c**2
    omega = sqrt( GM / (r * rschw)**3 )

    mdot_crit = 4 * pi * GM / (cgs_c * cgs_kapes)
    mdot_edd = 12 * mdot_crit
    facc = 3 * GM * (mdot * mdot_edd) / (8 * pi * (r * rschw)**3) &
    &    * (1 - sqrt(3/r))

    teff = ( facc / cgs_stef ) ** (1d0 / 4d0)
    zscale = sqrt( 2 * cgs_k_over_mh * teff ) / omega
  end subroutine

  elemental function fzscale(mbh, mdot, r) result(zscale)
    real(r64), intent(in) :: mbh, mdot, r
    real(r64) :: omega, facc, teff, zscale, rschw
    call cylinder(mbh, mdot, r, rschw, omega, facc, teff, zscale)
  end function

  elemental function fTeff(mbh, mdot, r) result(teff)
    real(r64), intent(in) :: mbh, mdot, r
    real(r64) :: omega, facc, teff, zscale, rschw
    call cylinder(mbh, mdot, r, rschw, omega, facc, teff, zscale)
  end function

  !--------------------------------------------------------------------------!

  elemental function kapabs0(X, Z) result(kab0)
    real(r64), intent(in) :: X, Z
    real(r64) :: kab0, kff0, kbf0
    kff0 = 3.68d22 * (1 - Z) * (1 + X)
    kbf0 = 4.34d25 * Z * (1 + X)
    kab0 = merge(kff0, 0.0_r64, use_opacity_ff)   &
    + merge(kbf0, 0.0_r64, use_opacity_bf)
  end function

  elemental function kapabp0(X, Z) result(kab0)
    real(r64), intent(in) :: X, Z
    real(r64) :: kab0
    kab0 = kapabs0(X, Z) * merge(37, 1, use_opacity_planck)
  end function

  elemental function kapes0(X, Z) result(kes0)
    real(r64), intent(in) :: X, Z
    real(r64) :: kes0
    kes0 = cgs_kapes_hydrogen * (1 + X) / 2 + 0*Z
  end function

  !--------------------------------------------------------------------------!

  elemental function fksct(rho,T) result(ksct)
    real(r64), intent(in) :: rho,T
    real(r64):: ksct0, ksct !, red
    ksct0 = kapes0(abuX,abuZ)
    ! red = 1d0 / (( 1 + 2.7e11 * rho / T**2 ) * ( 1 + (T / 4.5e8)**0.86 ))
    ksct = ksct0 ! * red
  end function

  elemental subroutine KAPPSCT(rho,T,ksct,krho,kT)
    real(r64), intent(in) :: rho,T
    real(r64), intent(out) :: ksct,krho,kT
    real(r64):: ksct0 !, red, redrho, redT
    ksct0 = kapes0(abuX,abuZ)
    ! red = 1d0 / (( 1 + 2.7e11 * rho / T**2 ) * ( 1 + (T / 4.5e8)**0.86 ))
    ksct = ksct0 ! * red
    !redrho = -2.7d+11/(T**2*(1.0d0 + 2.7d+11*rho/T**2)**2 &
    !    *(3.616073d-8*T**0.86d0 + 1.0d0))
    krho = 0 ! ksct0 * redrho
    !redT = -3.109823d-8*T**(-0.14d0)/((1.0d0 + 2.7d+11*rho/T**2)*( &
    !    3.616073d-8*T**0.86d0 + 1.0d0)**2) + 5.4d+11*rho/(T**3*(1.0d0 + &
    !    2.7d+11*rho/T**2)**2*(3.616073d-8*T**0.86d0 + 1.0d0))
    kT = 0 ! ksct0 * redT
  end subroutine

  !--------------------------------------------------------------------------!

  elemental function fkabs(rho,T) result(kabs)
    real(r64), intent(in) :: rho,T
    real(r64) :: kabs
    kabs = kapabs0(abuX,abuZ) * rho * T ** (-7d0/2d0)
  end function

  elemental subroutine kappabs(rho,T,kap,krho,kT)
    use ieee_arithmetic, only: ieee_is_nan
    real(r64), intent(in) :: rho,T
    real(r64), intent(out) :: kap,krho,kT

    kap = kapabs0(abuX,abuZ) * rho * T ** (-7d0/2d0)
    krho = kap / rho
    kT = (-7d0/2d0) * kap / T
  end subroutine

  !--------------------------------------------------------------------------!

  elemental function fkabp(rho,T) result(kabs)
    real(r64), intent(in) :: rho,T
    real(r64) :: kabs
    kabs = kapabp0(abuX,abuZ) * rho * T ** (-7d0/2d0)
  end function

  elemental subroutine kappabp(rho,T,kap,krho,kT)
    use ieee_arithmetic, only: ieee_is_nan
    real(r64), intent(in) :: rho,T
    real(r64), intent(out) :: kap,krho,kT

    kap = kapabp0(abuX,abuZ) * rho * T ** (-7d0/2d0)
    krho = kap / rho
    kT = (-7d0/2d0) * kap / T
  end subroutine

  !--------------------------------------------------------------------------!

  elemental function fkcnd(rho,T) result(kcnd)
    real(r64), intent(in) :: rho,T
    real(r64) :: kcnd
    kcnd = 5.6d-7 * cndredu * T**2.5d0
  end function

  elemental subroutine kappcnd(rho,T,kap,krho,kT)
    use ieee_arithmetic, only: ieee_is_nan
    real(r64), intent(in) :: rho,T
    real(r64), intent(out) :: kap,krho,kT

    kap = 5.6d-7 * cndredu * T**2.5d0
    krho = 0
    kT = 2.5d0 * kap / T
  end subroutine

end module
