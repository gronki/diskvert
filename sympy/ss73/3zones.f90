rho1 = 2.30485d-6*1d0/mbh*r**1.5d0/(alpha*cgs_kapes**3*mdot**2*f(r)**2)
T1 = 4.61893d+7*alpha**(-0.25d0)*cgs_kapes**(-0.25d0)*mbh**(-0.25d0)*r** &
      (-0.375d0)
H1 = 982981.0d0*cgs_kapes*mbh**1.0d0*mdot*f(r)
rho2 = 21.9208d0*alpha**(-0.7d0)*cgs_kapes**(-0.3d0)*mbh**(-0.7d0)*mdot &
      **0.4d0*r**(-1.65d0)*f(r)**0.4d0
T2 = 6.72321d+8*alpha**(-0.2d0)*cgs_kapes**0.2d0*mbh**(-0.2d0)*mdot** &
      0.4d0*r**(-0.9d0)*f(r)**0.4d0
H2 = 4639.53d0*alpha**(-0.1d0)*cgs_kapes**0.1d0*mbh**0.9d0*mdot**0.2d0*r &
      **1.05d0*f(r)**0.2d0
rho3 = 594575.0d0*alpha**(-0.7d0)*kappa_abs_0**(-0.15d0)*mbh**(-0.7d0)* &
      mdot**0.55d0*r**(-1.875d0)*f(r)**0.55d0
T3 = 744750.0d0*alpha**(-0.2d0)*kappa_abs_0**0.1d0*mbh**(-0.2d0)*mdot** &
      0.3d0*r**(-0.75d0)*f(r)**0.3d0
H3 = 154.415d0*alpha**(-0.1d0)*kappa_abs_0**0.05d0*mbh**0.9d0*mdot** &
      0.15d0*r**1.125d0*f(r)**0.15d0
r12 = 164.175d0*alpha**0.0952381d0*cgs_kapes**0.857143d0*mbh** &
      0.0952381d0*mdot**0.761905d0*f(r12)**0.761905d0
r23 = 5.05562d+19*cgs_kapes**1.33333d0*kappa_abs_0**(-0.666667d0)*mdot** &
      0.666667d0*f(r23)**0.666667d0
