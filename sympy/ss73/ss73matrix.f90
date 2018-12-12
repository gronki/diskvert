A(1, 1) = -4.0d0*H**3*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*rho*sqrt( &
      sol_mass)*sol_rschw**(-1.5d0)/mbh + mbh*mdot*sol_mdot_edd*f(r)
A(2, 1) = -0.5625d0*H**4*alpha*cgs_graw**1.5d0*r**(-4.5d0)*rho**2* &
      sol_mass**1.5d0*sol_rschw**(-4.5d0)*(T**(-3.5d0)*kappa_abs_0*rho &
      + kappa_es)/mbh**3 + 4.0d0*T**4*cgs_stef
A(3, 1) = H**2*cgs_graw*r**(-3.0d0)*rho*sol_mass/(mbh**2*sol_rschw**3) - &
      1.33333d0*T**4*cgs_stef/cgs_c - 2.0d0*T*cgs_boltz*rho/cgs_mhydr
M(1, 1) = 4.0d0*H**3*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*sqrt(sol_mass)* &
      sol_rschw**(-1.5d0)/mbh
M(2, 1) = 0.5625d0*H**4*T**(-3.5d0)*alpha*cgs_graw**1.5d0*kappa_abs_0*r &
      **(-4.5d0)*rho**2*sol_mass**1.5d0*sol_rschw**(-4.5d0)/mbh**3 + &
      1.125d0*H**4*alpha*cgs_graw**1.5d0*r**(-4.5d0)*rho*sol_mass** &
      1.5d0*sol_rschw**(-4.5d0)*(T**(-3.5d0)*kappa_abs_0*rho + kappa_es &
      )/mbh**3
M(3, 1) = -H**2*cgs_graw*r**(-3.0d0)*sol_mass/(mbh**2*sol_rschw**3) + &
      2.0d0*T*cgs_boltz/cgs_mhydr
M(1, 2) = 0
M(2, 2) = -1.96875d0*H**4*T**(-4.5d0)*alpha*cgs_graw**1.5d0*kappa_abs_0* &
      r**(-4.5d0)*rho**3*sol_mass**1.5d0*sol_rschw**(-4.5d0)/mbh**3 - &
      16.0d0*T**3*cgs_stef
M(3, 2) = 5.33333d0*T**3*cgs_stef/cgs_c + 2.0d0*cgs_boltz*rho/cgs_mhydr
M(1, 3) = 12.0d0*H**2*alpha*sqrt(cgs_graw)*pi*r**(-1.5d0)*rho*sqrt( &
      sol_mass)*sol_rschw**(-1.5d0)/mbh
M(2, 3) = 2.25d0*H**3*alpha*cgs_graw**1.5d0*r**(-4.5d0)*rho**2*sol_mass &
      **1.5d0*sol_rschw**(-4.5d0)*(T**(-3.5d0)*kappa_abs_0*rho + &
      kappa_es)/mbh**3
M(3, 3) = -2.0d0*H*cgs_graw*r**(-3.0d0)*rho*sol_mass/(mbh**2*sol_rschw** &
      3)
