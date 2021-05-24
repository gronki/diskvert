module globals

  use iso_fortran_env, only: r64 => real64
  use slf_cgs

  implicit none !--------------------------------------------------------!

  character(len = *), parameter :: version = '200708'

  real(r64), parameter, private :: X0 = 0.68d0 / cgs_kapes_hydrogen - 1
  ! real(r64), parameter, private :: X0 = 0.7381
  real(r64), parameter, private :: Z0 = 0.0134

  real(r64) :: abuX = X0
  real(r64) :: abuZ = Z0

  real(r64) :: kappa_es = cgs_kapes_hydrogen * (1 + X0) / 2
  real(r64) :: kappa_ff_0 = 3.68d22 * (1 - Z0) * (1 + X0)
  real(r64) :: kappa_bf_0 = 4.34d25 * Z0 * (1 + X0)
  real(r64) :: kappa_abs_0 = 3.68d22 * (1 - Z0) * (1 + X0)
  real(r64) :: opacities_kill = 1

  logical :: use_klein_nishina = .false.
  logical :: use_opacity_ff = .true.
  logical :: use_opacity_bf = .true.
  logical :: use_conduction = .false.
  logical :: use_relcompt = .true.

  real(r64) :: mbh, mdot, rschw
  real(r64), parameter :: accretion_efficiency = 0.083333_r64
  real(r64), parameter :: radius_isco = 3

  real(r64) :: cndredu = 1

  real(r64), parameter :: miu = 0.5_r64
  real(r64), parameter :: pi = 4*atan(real(1,r64))

  !   którego równania użyc?
  character, parameter :: EQUATION_DIFFUSION = 'D', &
  EQUATION_BALANCE = 'C', EQUATION_COMPTON = 'W'
  character :: cfg_temperature_method = EQUATION_DIFFUSION

contains !-----------------------------------------------------------------!

  elemental function mdot_edd(mbh)
    real(r64), intent(in) :: mbh
    real(r64) :: mdot_edd, GM

    GM = cgs_graw * (mbh * sol_mass)
    mdot_edd = 4 * pi * GM / (cgs_c * cgs_kapes * accretion_efficiency)
  end function

  elemental subroutine cylinder(mbh, mdot, r, rschw, omega, facc, teff, zscale)
    real(r64), intent(in) :: mbh, mdot, r
    real(r64), intent(out) :: rschw, omega, facc, teff, zscale
    real(r64) :: GM

    GM = cgs_graw * (mbh * sol_mass)

    rschw = 2 * GM / cgs_c**2
    omega = sqrt( GM / (r * rschw)**3 )

    facc = 3 * GM * (mdot * mdot_edd(mbh)) / (8 * pi * (r * rschw)**3) &
    &    * (1 - sqrt(radius_isco / r))

    teff = (facc / cgs_stef)**0.25_r64
    zscale = sqrt( 2 * cgs_k_over_mh * teff ) / omega
  end subroutine

  elemental function fzscale(mbh, mdot, r) result(zscale)
    real(r64), intent(in) :: mbh, mdot, r
    real(r64) :: omega, facc, teff, zscale, rschw
    call cylinder(mbh, mdot, r, rschw, omega, facc, teff, zscale)
  end function

  elemental function fzscale_rad(mbh, mdot, r) result(zscale)
    real(r64), intent(in) :: mbh, mdot, r
    real(r64) :: omega, facc, teff, zscale, rschw
    call cylinder(mbh, mdot, r, rschw, omega, facc, teff, zscale)
    zscale = cgs_kapes * facc / (omega**2 * cgs_c)
  end function

  elemental function fTeff(mbh, mdot, r) result(teff)
    real(r64), intent(in) :: mbh, mdot, r
    real(r64) :: omega, facc, teff, zscale, rschw
    call cylinder(mbh, mdot, r, rschw, omega, facc, teff, zscale)
  end function

  !--------------------------------------------------------------------------!

  elemental function kff0(X, Z)
    real(r64), intent(in) :: X, Z
    real(r64) :: kff0
    kff0 = 3.68e22_r64 * (1 - Z) * (1 + X)
  end function

  elemental function kbf0(X, Z)
    real(r64), intent(in) :: X, Z
    real(r64) :: kbf0
    kbf0 = 4.34e25_r64 * Z * (1 + X)
  end function

  elemental function kram0(X, Z) result(kab0)
    real(r64), intent(in) :: X, Z
    real(r64) :: kab0
    kab0 = merge(kff0(X, Z), 0.0_r64, use_opacity_ff)   &
         + merge(kbf0(X, Z), 0.0_r64, use_opacity_bf)
  end function

  elemental function kram0p(X, Z) result(kab0)
    real(r64), intent(in) :: X, Z
    real(r64) :: kab0
    ! for energy balance we use only free-free
    kab0 = merge(kff0(X, Z) * 37.0_r64, 0.0_r64, use_opacity_ff)
  end function

  !--------------------------------------------------------------------------!

  elemental function fkesp(rho,T) result(ksct)
    real(r64), intent(in) :: rho,T
    real(r64):: ksct

    associate (ksct0 => cgs_kapes_hydrogen * (1 + abuX) / 2, &
        kill => opacities_kill)
      if (use_klein_nishina) then
        ksct = ksct0 / ((1 + 2.7d+11*kill*rho/T**2)*(3.6161d-8*T**0.86d0*kill + 1))
      else
        ksct = ksct0
      end if
    end associate
  end function

  elemental subroutine kappesp(rho, T, ksct, krho, ktemp)
    real(r64), intent(in) :: rho, T
    real(r64), intent(out) :: ksct, krho, ktemp

    associate (ksct0 => cgs_kapes_hydrogen * (1 + abuX) / 2, &
        kill => opacities_kill)
      if (use_klein_nishina) then
        ksct = ksct0 / ((1 + 2.7d+11*kill*rho/T**2)*(3.6161d-8*T**0.86d0*kill + 1))
        krho = -2.7d+11*kill*ksct0/(T**2*(1.0d0 + 2.7d+11*kill*rho/T**2)**2*( &
              3.6161d-8*T**0.86d0*kill + 1.0d0))
        ktemp = -3.1098d-8*T**(-0.14d0)*kill*ksct0/((1.0d0 + 2.7d+11*kill*rho/T &
              **2)*(3.6161d-8*T**0.86d0*kill + 1.0d0)**2) + 5.4d+11*kill*ksct0* &
              rho/(T**3*(1.0d0 + 2.7d+11*kill*rho/T**2)**2*(3.6161d-8*T**0.86d0 &
              *kill + 1.0d0))
      else
        ksct = ksct0
        krho = 0
        ktemp = 0
      end if
    end associate
  end subroutine

  !--------------------------------------------------------------------------!

  elemental function fkabs(rho,T) result(kap)
    real(r64), intent(in) :: rho,T
    real(r64) :: kap

    associate (kbff0 => kram0(abuX,abuZ), kill => opacities_kill)
      kap = fkesp(rho, T) + kill * kbff0 * rho * T**(-3.5_r64)
    end associate
  end function

  elemental subroutine kappabs(rho,T,kap,krho,ktemp)
    use ieee_arithmetic, only: ieee_is_nan
    real(r64), intent(in) :: rho,T
    real(r64), intent(out) :: kap,krho,ktemp
    real(r64) :: kes, kes_rho, kes_temp

    call kappesp(rho, T, kes, kes_rho, kes_temp)

    associate (kbff0 => kram0(abuX,abuZ), kill => opacities_kill)
      kap = kes + kill * kbff0 * rho * T**(-3.5_r64)
      krho = kes_rho + kill * kbff0 * T**(-3.5_r64)
      ktemp = kes_temp - 3.5_r64 * kill * kbff0 * rho * T**(-4.5_r64)
    end associate
  end subroutine

  !--------------------------------------------------------------------------!

  elemental function fkabp(rho,T) result(kap)
    real(r64), intent(in) :: rho,T
    real(r64) :: kap

    associate (kbff0 => kram0p(abuX,abuZ))
      kap = kbff0 * rho * T**(-3.5_r64)
    end associate
  end function

  elemental subroutine kappabp(rho,T,kap,krho,ktemp)
    use ieee_arithmetic, only: ieee_is_nan
    real(r64), intent(in) :: rho,T
    real(r64), intent(out) :: kap,krho,ktemp
    
    associate (kbff0 => kram0p(abuX,abuZ))
      kap = kbff0 * rho * T**(-3.5_r64)
      krho = kbff0 * T**(-3.5_r64)
      ktemp = - 3.5_r64 * kbff0 * rho * T**(-4.5_r64)
    end associate
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

  !--------------------------------------------------------------------------!

  elemental function rho_crit(Trad, T) result(rho)
    real(r64), intent(in) :: Trad
    real(r64), intent(in), optional :: T
    real(r64) :: rho
    
    rho = cgs_kapes / kram0p(abuX, abuZ) &
          * 4 * cgs_k_over_mec2 * Trad**4.5
    if (present(T)) then
      rho = rho * sqrt(Trad / T)
    end if
  end function

end module
