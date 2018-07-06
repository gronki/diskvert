module ss73solution

    use iso_fortran_env, only: r64 => real64
    use slf_cgs
    use globals, only: kappa_es, kappa_abs_0, pi

    implicit none

    private :: f

contains

!--------------------------------- APXDISK2D ----------------------------------!
!      calculates a radial slice through an accretion disk, based on SS73      !
!                             analytical solution                              !
!----------------------------------- INPUTS -----------------------------------!
!             mbh: black hole mass                                             !
!            mdot: accretion rate                                              !
!               r: array of radii                                              !
!           alpha: alpha parameter                                             !
!               z: array of heights to calculate, expressed in rschw           !
!---------------------------------- OUTPUTS -----------------------------------!
!                               rho: density                                   !
!                                 T: temperature                               !
!------------------------------------------------------------------------------!

  subroutine apxdisk2d(mbh, mdot, r, alpha, z, rho, T, Frad)
    use globals, only: fteff
    real(r64), intent(in) :: mbh, mdot,  alpha
    real(r64), intent(in), dimension(:) :: r, z
    real(r64), intent(out), dimension(:,:) :: rho,T,Frad
    real(r64) :: r12, r23, rschw
    real(r64), dimension(size(r)) :: rhoc, Tc, Teff, H
    integer :: i,j

    teff = fteff(mbh,mdot,r)
    rschw = sol_rschw * mbh

    call ss73_zonebounds(mbh, mdot, alpha, r12, r23)
    call ss73_sel(mbh, mdot, r, alpha, r12, r23, rhoc, Tc, H)

    do i = 1,size(R)
      call ss73_refin(mbh, mdot, r(i), alpha, rhoc(i), Tc(i), H(i))
    end do

    do concurrent (i = 1:size(z), j = 1:size(r))
      call apxdisk(rhoc(j), Tc(j), Teff(j), H(j), rschw*z(i),  &
            & rho(i,j), T(i,j), Frad(i,j))
    end do
  end subroutine

  subroutine apxdisk2d_c(mbh, mdot, r, alpha, z, rho, T, Frad, nr, nz)  &
        &  bind(C, name = 'apxdisk2d')
    use iso_c_binding, only: c_double, c_int
    integer(c_int), value :: nr,nz
    real(c_double), intent(in), value :: mbh, mdot,  alpha
    real(c_double), intent(in) :: r(nr), z(nz)
    real(c_double), intent(out) :: rho(nz,nr), T(nz,nr), Frad(nz,nr)

    call apxdisk2d(mbh, mdot, r, alpha, z, rho, T, Frad)
  end subroutine

!---------------------------------- APXDISK -----------------------------------!
!  for given parameters (obtained by another routines) helps to calculate the  !
!    approximate vertical structure of an accretion disk, assuming gaussian    !
!                 vertical profile of density and temperature                  !
!----------------------------------- INPUTS -----------------------------------!
!                 rhoc: central density                                        !
!                   Tc: central temperature                                    !
!                 Teff: effective temperature                                  !
!                    H: disk characteristic height                             !
!                    z: height to calculate parameters (in cm)                 !
!---------------------------------- OUTPUTS -----------------------------------!
!                               rho: density                                   !
!                                 T: temperature                               !
!------------------------------------------------------------------------------!

  elemental subroutine apxdisk(rhoc, Tc, Teff, H, z, rho, T, Frad)
    use slf_cgs, only: cgs_stef
    real(r64), intent(in) :: rhoc, Tc, Teff, H, z
    real(r64), intent(out) :: rho, T, Frad
    real(r64) :: expo, x
    expo = exp(- (z/H)**2 / 2)
    rho = rhoc * expo
    T = (Tc - Teff * 0.841) * expo + Teff * 0.841
    x = merge(0.5*z/H, 1.0_r64, z < 2*H)
    Frad = (2 - X) * X * cgs_stef * Teff**4
  end subroutine

!------------------------------------- F --------------------------------------!
!                       computes the 1-sqrt(3/r) factor                        !
!----------------------------------- INPUTS -----------------------------------!
!                          r: radius (dimensionless)                           !
!---------------------------------- RETURNS -----------------------------------!
!                      f-factor for r > 3, zero otherwise                      !
!------------------------------------------------------------------------------!

    elemental function f(r)
        real(r64), intent(in) :: r
        real(r64) :: f
        if ( r > 3 ) then
            f = 1 - sqrt(3/r)
        else
            f = 0
        end if
    end function

!------------------------------ DVss73_ZONEBOUNDS ------------------------------!
!   computes the boundaries between Prad+elsct, Pgas+elsc, Pgas+kapff zones    !
!---------------------------------- OUTPUTS -----------------------------------!
!                  r12: boundary radius between zones 1 and 2                  !
!                  r23: boundary radius between zones 2 and 3                  !
!------------------------------------------------------------------------------!

    elemental subroutine ss73_zonebounds(mbh, mdot, alpha, r12,r23)
        real(r64), intent(in) :: mbh, mdot, alpha
        real(r64), intent(out) :: r12,r23
        integer :: i
        r12 = 164.17d0 * alpha**0.095238d0 * kappa_es**0.85714d0 &
                & * mbh**0.095238d0 * mdot**0.7619d0
        r23 = 5.0556d+19 * kappa_es**1.3333d0 * kappa_abs_0**(-0.66667d0) &
                & * mdot**0.66667d0
        do i = 1,8
            r12 = 164.17d0 * alpha**0.095238d0 * kappa_es**0.85714d0 &
                & * mbh**0.095238d0 * mdot**0.7619d0 * f(r12)**0.7619d0
            r23 = 5.0556d+19 * kappa_es**1.3333d0 * kappa_abs_0**(-0.66667d0) &
                & * mdot**0.66667d0 * f(r23)**0.66667d0
        end do
    end subroutine

!-------------------------------- DVss73_ZONE1 ---------------------------------!
!     calculates rho, T, H for radiation pressure and electron scattering      !
!                              dominated zone (1)                              !
!----------------------------------- INPUTS -----------------------------------!
!                                  r: radius                                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    elemental subroutine ss73_zone1(mbh, mdot, r, alpha, rho, T, H)
        real(r64), intent(in) :: mbh, mdot, r, alpha
        real(r64), intent(out) :: rho,T,H
        rho = 2.3049d-6*1.0/mbh*r**1.5d0/(alpha*kappa_es**3*mdot**2*f(r)**2)
        T = 4.6189d+7*alpha**(-0.25d0)*kappa_es**(-0.25d0)*mbh**(-0.25d0)*r**( &
              -0.375d0)
        H = 9.8298d+5*kappa_es*mbh**1.0d0*mdot*f(r)
    end subroutine

!-------------------------------- DVss73_ZONE2 ---------------------------------!
! calculates rho, T, H for gas pressure and electron scattering dominated zone !
!                                     (2)                                      !
!----------------------------------- INPUTS -----------------------------------!
!                                  r: radius                                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    elemental subroutine ss73_zone2(mbh, mdot, r, alpha, rho, T, H)
        real(r64), intent(in) :: mbh, mdot, r, alpha
        real(r64), intent(out) :: rho,T,H
        rho = 21.921d0*alpha**(-0.7d0)*kappa_es**(-0.3d0)*mbh**(-0.7d0)*mdot** &
              0.4d0*r**(-1.65d0)*f(r)**0.4d0
        T = 6.7232d+8*alpha**(-0.2d0)*kappa_es**0.2d0*mbh**(-0.2d0)*mdot** &
              0.4d0*r**(-0.9d0)*f(r)**0.4d0
        H = 4639.5d0*alpha**(-0.1d0)*kappa_es**0.1d0*mbh**0.9d0*mdot**0.2d0*r &
              **1.05d0*f(r)**0.2d0
    end subroutine

!-------------------------------- DVAXP_ZONE3 ---------------------------------!
!  calculates rho, T, H for gas pressure and free-free opacity dominated zone  !
!                                     (3)                                      !
!----------------------------------- INPUTS -----------------------------------!
!                                  r: radius                                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    elemental subroutine ss73_zone3(mbh, mdot, r, alpha, rho, T, H)
        real(r64), intent(in) :: mbh, mdot, r, alpha
        real(r64), intent(out) :: rho,T,H
        rho = 5.9458d+5*alpha**(-0.7d0)*kappa_abs_0**(-0.15d0)*mbh**(-0.7d0)*mdot &
              **0.55d0*r**(-1.875d0)*f(r)**0.55d0
        T = 7.4475d+5*alpha**(-0.2d0)*kappa_abs_0**0.1d0*mbh**(-0.2d0)*mdot**0.3d0 &
              *r**(-0.75d0)*f(r)**0.3d0
        H = 154.42d0*alpha**(-0.1d0)*kappa_abs_0**0.05d0*mbh**0.9d0*mdot**0.15d0*r &
              **1.125d0*f(r)**0.15d0
    end subroutine

!---------------------------------- ss73_SEL -----------------------------------!
!     selects the zone based on given radius and computes appropriate disk     !
!                                  parameters                                  !
!----------------------------------- INPUTS -----------------------------------!
!                           r: radius                                          !
!                         r12: radius of 1-2 boundary                          !
!                         r23: radius of 2-3 boundary                          !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    elemental subroutine ss73_sel(mbh, mdot, r, alpha, r12, r23, rho, T, H)
        real(r64), intent(in) :: mbh, mdot, r, alpha, r12, r23
        real(r64), intent(out) :: rho, T, H

        if (r < r12) then
            call ss73_zone1(mbh, mdot, r, alpha, rho, T, H)
        else if (r > r23) then
            call ss73_zone3(mbh, mdot, r, alpha, rho, T, H)
        else
            call ss73_zone2(mbh, mdot, r, alpha, rho, T, H)
        end if
    end subroutine

!--------------------------------- ss73_ESTIM ----------------------------------!
! Computes disk density, temperature and height based on 3-zone approximation. !
!        Calls DVss73_ZONEBOUNDS to calculate boundaries between zones.         !
!----------------------------------- INPUTS -----------------------------------!
!                                  r: radius                                   !
!---------------------------------- OUTPUTS -----------------------------------!
!                           rho: central density                               !
!                             T: central temperature                           !
!                             H: disk height                                   !
!------------------------------------------------------------------------------!

    pure subroutine ss73_estim(mbh, mdot, r, alpha, rho, T, H)
        real(r64), intent(in) :: mbh, mdot, r, alpha
        real(r64), intent(out) :: rho, T, H
        real(r64) :: r12, r23

        call ss73_zonebounds(mbh, mdot, alpha, r12, r23)
        call ss73_sel(mbh, mdot, r, alpha, r12, r23, rho, T, H)

    end subroutine

!--------------------------------- ss73_MATRIX ---------------------------------!
!     calculates matrix and RHS of the equation system for solving general     !
!                           accretion disk structure                           !
!----------------------------------- INPUTS -----------------------------------!
!                                 r: radius                                    !
!                               rho: density                                   !
!                                 T: temperature                               !
!                                 H: disk height                               !
!---------------------------------- OUTPUTS -----------------------------------!
!                           A: RHS of the equations                            !
!                           M: coefficient matrix                              !
!------------------------------------------------------------------------------!

  pure subroutine ss73_matrix(mbh, mdot, r, alpha, rho, T, H, A, M)
    real(r64), intent(in) :: mbh, mdot, r, alpha, rho, T, H
    real(r64), intent(out) :: A(3,1), M(3,3)

    A(1, 1) = -4*H**3*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*rho*sqrt(sol_mass) &
          /(mbh*sol_rschw**(3.0d0/2.0d0)) + mbh*mdot*sol_mdot_edd*( &
          -1.73205080756888d0*r**(-0.5d0) + 1)
    A(2, 1) = -9.0d0/16.0d0*H**4*alpha*cgs_graw**(3.0d0/2.0d0)*r**(-4.5d0)* &
          rho**2*sol_mass**(3.0d0/2.0d0)*(kappa_es + kappa_abs_0*rho/T**(7.0d0/ &
          2.0d0))/(mbh**3*sol_rschw**(9.0d0/2.0d0)) + 4*T**4*cgs_stef
    A(3, 1) = H**2*cgs_graw*r**(-3.0d0)*rho*sol_mass/(mbh**2*sol_rschw**3) &
          - 4.0d0/3.0d0*T**4*cgs_stef/cgs_c - 2*T*cgs_boltz*rho/cgs_mhydr
    M(1, 1) = 4*H**3*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*sqrt(sol_mass)/( &
          mbh*sol_rschw**(3.0d0/2.0d0))
    M(2, 1) = (9.0d0/8.0d0)*H**4*alpha*cgs_graw**(3.0d0/2.0d0)*r**(-4.5d0)* &
          rho*sol_mass**(3.0d0/2.0d0)*(kappa_es + kappa_abs_0*rho/T**(7.0d0/ &
          2.0d0))/(mbh**3*sol_rschw**(9.0d0/2.0d0)) + (9.0d0/16.0d0)*H**4* &
          alpha*cgs_graw**(3.0d0/2.0d0)*kappa_abs_0*r**(-4.5d0)*rho**2* &
          sol_mass**(3.0d0/2.0d0)/(T**(7.0d0/2.0d0)*mbh**3*sol_rschw**( &
          9.0d0/2.0d0))
    M(3, 1) = -H**2*cgs_graw*r**(-3.0d0)*sol_mass/(mbh**2*sol_rschw**3) + 2 &
          *T*cgs_boltz/cgs_mhydr
    M(1, 2) = 0
    M(2, 2) = -63.0d0/32.0d0*H**4*alpha*cgs_graw**(3.0d0/2.0d0)*kappa_abs_0*r** &
          (-4.5d0)*rho**3*sol_mass**(3.0d0/2.0d0)/(T**(9.0d0/2.0d0)*mbh**3 &
          *sol_rschw**(9.0d0/2.0d0)) - 16*T**3*cgs_stef
    M(3, 2) = (16.0d0/3.0d0)*T**3*cgs_stef/cgs_c + 2*cgs_boltz*rho/cgs_mhydr
    M(1, 3) = 12*H**2*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*rho*sqrt(sol_mass) &
          /(mbh*sol_rschw**(3.0d0/2.0d0))
    M(2, 3) = (9.0d0/4.0d0)*H**3*alpha*cgs_graw**(3.0d0/2.0d0)*r**(-4.5d0)* &
          rho**2*sol_mass**(3.0d0/2.0d0)*(kappa_es + kappa_abs_0*rho/T**(7.0d0/ &
          2.0d0))/(mbh**3*sol_rschw**(9.0d0/2.0d0))
    M(3, 3) = -2*H*cgs_graw*r**(-3.0d0)*rho*sol_mass/(mbh**2*sol_rschw**3)

  end subroutine

!--------------------------------- ss73_refin ---------------------------------!
!  solves the equation system for disk structure, starting from given initial  !
!                                   solution                                   !
!----------------------------------- INPUTS -----------------------------------!
!                             r: radius                                        !
!                           rho: initial density                               !
!                             T: initial temperature                           !
!                             H: initial disk height                           !
!---------------------------------- OUTPUTS -----------------------------------!
!                        rho: final central density                            !
!                          T: final central temperature                        !
!                          H: final disk height                                !
!------------------------------------------------------------------------------!

    subroutine ss73_refin(mbh, mdot, r, alpha, rho, T, H)
        real(r64), intent(in) :: mbh, mdot, r, alpha
        real(r64), intent(inout) :: rho,T,H
        real(r64), dimension(3,3) :: M
        real(r64), dimension(3) :: A
        integer, dimension(3) :: ipiv
        integer :: errno
        integer :: i
        integer, parameter :: niter = 24

        convergence: do i = 1,niter
            call ss73_matrix(mbh, mdot, r, alpha, rho, T, H, A, M)
            call dgesv(3, 1, M, 3, ipiv, A, 3, errno)
            if ( errno .ne. 0 ) then
                error stop "error while converging to the solution"
            end if
            rho = rho + A(1) * ramp(i,niter)
            T = T + A(2) * ramp(i,niter)
            H = H + A(3) * ramp(i,niter)
        end do convergence

    contains

        !----------------------------------------------------!
        !  This function generates a smooth, S-shaped curve  !
        !----------------------------------------------------!
        elemental function ramp(i,n) result(y)
            integer, intent(in) :: i,n
            real(r64) :: t,y
            t = merge(real(i) / n, 1.0, i .le. n)
            y = 3 * t**2 - 2 * t**3
        end function

    end subroutine

end module
