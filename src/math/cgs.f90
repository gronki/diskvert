module slf_cgs

  use iso_fortran_env, only: r64 => real64
  implicit none

  real(r64), parameter, private :: pi = 4*atan(real(1,r64))

  real(r64), parameter :: &
    cgs_boltz = 1.3806581212d-16,     &
    cgs_c = 2.99792458d10, &
    cgs_h = 6.62607554040d-27, &
    cgs_graw = 6.67259858585d-8, &
    cgs_mel = 9.10938975454d-28,    &
    cgs_hbar= cgs_h / (2*pi), &
    cgs_stef = 2 * pi**5 * cgs_boltz**4 / ( 15 * cgs_h**3 * cgs_c**2 ), &
    cgs_a = 4 * cgs_stef / cgs_c, &
    cgs_alpha = 1 / 137.036d0, &
    cgs_qel = sqrt( cgs_alpha * cgs_hbar * cgs_c ), &
    cgs_thomson = 8 * pi / 3 * cgs_qel**4 / ( cgs_c**4 * cgs_mel**2 ), &
    cgs_elrad = cgs_qel**2 / ( cgs_c**2 * cgs_mel ), &
    cgs_mhydr = 1.6733d-24,    &
    cgs_kapes_hydrogen = cgs_thomson / cgs_mhydr, &
    cgs_k = cgs_boltz,   &
    cgs_k_over_mh = cgs_boltz / cgs_mhydr , &
    cgs_mh_over_k = cgs_mhydr / cgs_boltz, &
    cgs_k_over_mec = cgs_boltz / ( cgs_mel * cgs_c ), &
    cgs_mec_over_k = ( cgs_mel * cgs_c ) / cgs_boltz, &
    cgs_k_over_mec2 = cgs_boltz / ( cgs_mel * cgs_c**2 ), &
    cgs_mec2_over_k = ( cgs_mel * cgs_c**2 ) / cgs_boltz, &
    cgs_msun = 1.98855d33, &
    cgs_lsun = 3.828d33,     &
    cgs_kapes = cgs_kapes_hydrogen * 0.85

  real(r64), parameter :: &
    keV_in_erg = 1.6021772d-9, &
    keV_in_kelvin = cgs_boltz / keV_in_erg,   &
    angstr_keV = 1e8 * cgs_h * cgs_c / keV_in_erg

  ! UWAGA: nie uzywac!!!!
  real(r64), parameter ::  &
    sol_mass = 1.98855d33, sol_lum = 3.828d33,              &
    sol_rschw = 2 * cgs_graw * sol_mass / cgs_c**2,         &
    sol_mdot_crit = 4 * pi * (cgs_graw * sol_mass) / (cgs_c * cgs_kapes), &
    sol_mdot_edd = 12 * sol_mdot_crit,                      &
    sol_facc_edd = 3 * cgs_graw * sol_mass * sol_mdot_edd   &
              / ( 8 * pi * sol_rschw**3 ),            &
    sol_omega = sqrt( cgs_graw * sol_mass / sol_rschw**3 )

end module
