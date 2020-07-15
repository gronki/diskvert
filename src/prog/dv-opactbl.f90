program dvopa

  use globals
  use slf_cgs
  use settings, only: rdargvgl, outfn
  use ranges, only: logrange
  use heatbalance, only: fcool2

  implicit none
  real(r64) :: trad, wsp1, wsp2, wsp3
  integer :: un, i, j
  integer, parameter :: n = 256
  real(r64) :: rho(n), T(n) , rho0, cool(n), cool_dr(n), cool_dT(n), cool_p0(n)

  outfn = ""
  call rdargvgl

  write(*, '(a)', advance = 'no') "input T_rad = "
  read(*, *) trad

  write(*, *)
  
  wsp1 = (8. / 3.) * cgs_kapes / kram0p(abuX, abuZ) * cgs_k_over_mec2 * 1e27_r64
  wsp2 = 4 * cgs_stef * (cgs_kapes * 4 * cgs_k_over_mec2)**2 / kram0p(abuX, abuZ) * 1e57_r64
  write (*, '(a15, 2es9.2)') 'coeff = ', wsp1, wsp1 / cgs_mhydr
  write (*, '(a15, es9.2)') 'coeff2 = ', wsp2
  write (*, '(a15, es9.2)') 'cool0 = ', wsp2 * (Trad / 1e6)**9.5 

  rho0 = wsp1 * (trad / 1e6)**4.5
  write (*, '(a15, es9.2)') 'rho0 = ', rho0 / cgs_mhydr
  write (*, '(a15, es9.2)') 'kram0 = ', kram0(abuX, abuZ)
  write (*, '(a15, es9.2)') 'kram0p = ', kram0p(abuX, abuZ)

  call logrange(T, Trad, Trad * 3e4)
  call logrange(rho, rho0 * 3e-5, rho0 * 3e4)

  if (outfn == "") write(outfn, '("cool_logT_", f0.2, ".dat")') log10(trad)
  open(newunit = un, file = outfn, action = 'write')

  write(un, '(i8, 3es12.4)') n, trad, rho0 / cgs_mhydr, wsp2 * (Trad / 1e6)**9.5 
  write(un, *)
  write(un, '(*(es12.4))') T(:)
  write(un, '(*(es12.4))') rho(:) / cgs_mhydr

  do j = 1, n
    call fcool2(rho(j), T, trad, cool, cool_dr, cool_dT)
    cool_dr(:) = cool_dr * rho(j) / cool
    cool_dT(:) = cool_dT * T      / cool
    cool_p0(:) = cool_dT - cool_dr
    ! cool = cool / wsp2
    write(un, '("# j =", i4, ", rho =", es12.4)') j, rho(j) / cgs_mhydr
    write(un, '(*(es12.4))') cool
    write(un, '(*(es12.4))') cool_dr
    write(un, '(*(es12.4))') cool_dT
    write(un, '(*(es12.4))') cool_p0
  end do

  close(un)

  write(*, '(a, a)') 'output: ', trim(outfn)

end program