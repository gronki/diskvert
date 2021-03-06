module heatbalance


  use slf_cgs
  use globals, only: kappabp, kappesp, miu
  use slf_threshold ! <--- TODO: fix fortdep
  use iso_fortran_env, only: r64 => real64
  implicit none

  private :: ramp

contains

  !----------------------------------------------------------------------------!
  ! THIS IS DEPRECATED! do not use! use heatbil2 instead

  pure subroutine heatbil(T, Trad, Pgas, heat)

    real(r64), intent(inout) :: T
    real(r64), intent(in) :: Trad, Pgas, heat
    real(r64), dimension(3,2) :: kap
    real(r64), dimension(2,2) :: F
    real(r64) :: rho, bil, bil_T, rhotemp
    integer :: i
    integer, parameter :: niter = 64

    rhotemp = Pgas * miu * cgs_mhydr / cgs_boltz

    converge: do i = 1,niter

      rho = rhotemp / T

      call kappabp(rho, T, kap(1,1), kap(2,1), kap(3,1))
      call kappesp(rho, T, kap(1,2), kap(2,2), kap(3,2))

      F(1,:) = kap(1,:)
      F(2,:) = kap(3,:) - rho / T * kap(2,:)

      bil = (-T*heat + 4*cgs_stef*rhotemp*(T - Trad)*(4*Trad**4* &
      cgs_k_over_mec2*F(1, 2) + (T + Trad)*(T**2 + Trad**2)*F(1, 1)))/( &
      T*Trad**4)
      bil_T = 4*cgs_stef*rhotemp*(T**5*F(2, 1) + 3*T**4*F(1, 1) + 4*T**2*Trad &
            **4*cgs_k_over_mec2*F(2, 2) - 4*T*Trad**5*cgs_k_over_mec2*F(2, 2 &
            ) - T*Trad**4*F(2, 1) + 4*Trad**5*cgs_k_over_mec2*F(1, 2) + Trad &
            **4*F(1, 1))/(T**2*Trad**4)


      T = T - ramp(i,niter) * bil / bil_T

    end do converge

  end subroutine

  !----------------------------------------------------------------------------!
  ! cooling funciton for heatbil2

  elemental subroutine fcool2(rho, temp, trad, cool, cool_drho, cool_dtemp)
    use globals, only: use_relcompt
    real(r64), intent(in) :: rho, temp, trad
    real(r64), intent(out) :: cool
    real(r64), intent(out), optional :: cool_drho, cool_dtemp
    real(r64) :: kabpv(3), ksctv(3)

    call kappesp(rho, temp, ksctv(1), ksctv(2), ksctv(3))
    call kappabp(rho, temp, kabpv(1), kabpv(2), kabpv(3))

    cool = 4*cgs_stef*rho*(4*cgs_k_over_mec2*trad**4*(temp*merge(4* &
      cgs_k_over_mec2*temp + 1, 1.0d0, use_relcompt) - trad)* &
      ksctv(1) + (temp**4 - trad**4)*kabpv(1))

    if (present(cool_drho)) then
      if (use_relcompt) then
        cool_drho = 4*cgs_stef*(4*cgs_k_over_mec2*trad**4*(temp*(4* &
          cgs_k_over_mec2*temp + 1) - trad)*ksctv(1) + rho*(4* &
          cgs_k_over_mec2*trad**4*(temp*(4*cgs_k_over_mec2*temp + 1) - trad &
          )*ksctv(2) + (temp**4 - trad**4)*kabpv(2)) + (temp**4 - trad**4)* &
          kabpv(1))
      else
        cool_drho = 4*cgs_stef*(4*cgs_k_over_mec2*trad**4*(1.0d0*temp - trad) &
          *ksctv(1) + rho*(4*cgs_k_over_mec2*trad**4*(1.0d0*temp - trad)* &
          ksctv(2) + (temp**4 - trad**4)*kabpv(2)) + (temp**4 - trad**4)* &
          kabpv(1))
      end if
    end if

    if (present(cool_dtemp)) then
      if (use_relcompt) then
        cool_dtemp = 4*cgs_stef*rho*(4*cgs_k_over_mec2*trad**4*(8* &
          cgs_k_over_mec2*temp + 1)*ksctv(1) + 4*cgs_k_over_mec2*trad**4*( &
          temp*(4*cgs_k_over_mec2*temp + 1) - trad)*ksctv(3) + 4*temp**3* &
          kabpv(1) + (temp**4 - trad**4)*kabpv(3))
      else
        cool_dtemp = 4*cgs_stef*rho*(4*cgs_k_over_mec2*trad**4*(1.0d0*temp - &
          trad)*ksctv(3) + 4.0d0*cgs_k_over_mec2*trad**4*ksctv(1) + 4*temp &
          **3*kabpv(1) + (temp**4 - trad**4)*kabpv(3))
      end if
    end if
  end subroutine

  !----------------------------------------------------------------------------!
  ! solves the heating-cooling balance equation

  pure subroutine heatbil2(rho, temp, trad, heat, isobar)

    real(r64), intent(in) :: trad, heat
    real(r64), intent(inout) :: rho, temp
    logical, intent(in) :: isobar
    real(r64) :: rhotemp, cool, cool_drho, cool_dtemp, dcool
    integer :: i
    integer, parameter :: niter = 36

    ! keep rho*temp value for isobaric case
    rhotemp = rho * temp

    ! some initial guess
    temp = heat  &
          & / (16 * cgs_k_over_mec2 * cgs_kapes * rho &
          & * cgs_stef * trad**4 )
    temp = sqrt(trad**2 + temp**2)

    do i = 1, niter
      if (.not. isobar) then
        ! for rho = const. phase, we only move temperature
        call fcool2(rho, temp, trad, cool, cool_dtemp = dcool)
      else
        ! for isobaric, we need both derivatives
        call fcool2(rho, temp, trad, cool, cool_drho = cool_drho, cool_dtemp = cool_dtemp)
        dcool = cool_dtemp - rho / temp * cool_drho
      end if

      ! apply the correction and change the density for isobaric case
      temp = temp - (cool - heat) / dcool * ramp(i,niter)
      if (isobar) rho = rhotemp / temp
    end do
  end subroutine heatbil2

  !----------------------------------------------------------------------------!

  elemental function ramp(i,n) result(y)
    integer, intent(in) :: i,n
    real(r64) :: t,y
    t = merge(real(i) / n, 1.0, i .le. n)
    y = (3*t**2 - 8*t + 6) * t**2
  end function

  !----------------------------------------------------------------------------!

end module


module energy_balance

    use iso_fortran_env, only: r64 => real64
    use globals
    use slf_threshold

    implicit none

    real(r64) ::   fcool_heat, &
                & fcool_pgas, &
                & fcool_prad, &
                & fcool_trad

contains


    pure subroutine fcool(logTemp,y,dy)

        real(r64), intent(in) :: logTemp
        real(r64), intent(out) :: y
        real(r64), intent(out), optional :: dy
        real(r64) :: xxA2kesPradk_over_mec, xx4sgT4_m_3cPrad, xx4sgT4, xxA
        real(r64) :: rho,kabs,cool
        real(r64) :: Temp

        Temp = exp(logTemp)

        rho = fcool_pgas * miu * cgs_mhydr / ( cgs_boltz * Temp )
        kabs = fkabs(rho,temp)

        xx4sgT4 = 4 * cgs_stef * Temp**4
        xx4sgT4_m_3cPrad = xx4sgT4 - 3 * cgs_c * fcool_Prad
        xxA2kesPradk_over_mec = 12 * kappa_es * fcool_prad * cgs_k_over_mec

        xxA = kabs * xx4sgT4_m_3cPrad + xxA2kesPradk_over_mec * ( Temp - fcool_Trad )

        cool = xxA * rho

        y = fcool_heat - cool

        if ( present(dy) ) then
            xxA = 4 * (xx4sgT4 - 9 * xx4sgT4_m_3cPrad / 8) * kabs &
                & + xxA2kesPradk_over_mec * Temp
            dy =  - ( rho * xxA - cool )
        end if

    end subroutine

    pure subroutine calc_compswitch(pgas, tgas, compsw)
        real(r64), intent(in) :: tgas,pgas
        real(r64), intent(out) ::compsw
        real(r64) :: rho, epsi, taues, compy, kabs

        rho =  pgas * miu * cgs_mhydr / ( cgs_boltz * tgas )

        kabs = fkabs(rho, tgas)
        epsi = kabs / (kabs + kappa_es)
        taues = (1 - epsi) / sqrt(epsi)
        compy = 4 * cgs_boltz * tgas &
            & / ( cgs_mel * cgs_c**2 ) &
            & * max( taues, taues**2 )
        compsw = thrtanh( log(compy) / log(2d0) )

    end subroutine

    pure subroutine calc_compW(heat, pgas, prad, heat_max, compw)
        real(r64), intent(in) :: heat, pgas, prad
        real(r64), intent(out) :: heat_max, compw

        heat_max = 12 * cgs_kapes * prad  &
            & * ( pgas * miu * cgs_mhydr ) &
            & / ( cgs_mel * cgs_c )

        compw = heat / ( heat_max - heat )
    end subroutine

end module
