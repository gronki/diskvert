cool = 4*cgs_stef*rho*(4*cgs_k_over_mec2*trad**4*(temp*merge(4* &
      cgs_k_over_mec2*temp + 1, 1.0d0, use_precise_balance) - trad)* &
      ksctv(1) + (temp**4 - trad**4)*kabpv(1))
if (use_precise_balance) then
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
if (use_precise_balance) then
   cool_dtemp = 4*cgs_stef*rho*(4*cgs_k_over_mec2*trad**4*(8* &
      cgs_k_over_mec2*temp + 1)*ksctv(1) + 4*cgs_k_over_mec2*trad**4*( &
      temp*(4*cgs_k_over_mec2*temp + 1) - trad)*ksctv(3) + 4*temp**3* &
      kabpv(1) + (temp**4 - trad**4)*kabpv(3))
else
   cool_dtemp = 4*cgs_stef*rho*(4*cgs_k_over_mec2*trad**4*(1.0d0*temp - &
      trad)*ksctv(3) + 4.0d0*cgs_k_over_mec2*trad**4*ksctv(1) + 4*temp &
      **3*kabpv(1) + (temp**4 - trad**4)*kabpv(3))
end if
