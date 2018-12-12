from sympy import fcode, symbols, Rational, Eq, var, Symbol, Matrix, \
    MatrixSymbol, sqrt, Function

r, rho, T, H = symbols('r, rho, T, H', real = True, positive = True)
kappa_es, kappa_ff = symbols('kappa_es kappa_abs_0', real = True, positive = True)
r, mbh, alpha, mdot, omega = symbols('r mbh alpha mdot omega', real = True, positive = True)
pi,cgs_graw,sol_mass,cgs_stef,cgs_boltz = symbols('pi cgs_graw sol_mass cgs_stef cgs_boltz', real = True, positive = True)
sol_rschw,sol_mdot_edd,cgs_mhydr,cgs_c = symbols('sol_rschw sol_mdot_edd cgs_mhydr cgs_c', real = True, positive = True)

mdot_cgs = mdot * mbh * sol_mdot_edd
omega = sqrt(cgs_graw * sol_mass / sol_rschw**3) / (mbh * r**1.5)

f = Function('f')

kappa = lambda rho,T: kappa_es + kappa_ff * rho * T**Rational(-7,2)
P = lambda rho,T: 2 * cgs_boltz * rho * T / cgs_mhydr \
    + 4 * cgs_stef / (3 * cgs_c) * T**4

eq1 = 4 * pi * alpha * rho * H**3 * omega - mdot_cgs * f(r)
eq2 = ( Rational(3,4) * rho * H**2 )**2 * alpha * omega**3 * kappa(rho,T) \
    - 4 * cgs_stef * T**4
eq3 = P(rho,T) - rho * H**2 * omega**2

A = Matrix([eq1,eq2,eq3])

def f90code(r,l):
    f90sym, f90fun, f90cod = fcode(r.evalf(6), l, source_format = 'free', \
        standard = 2008, human = False)
    return f90cod

with open('ss73matrix.f90','w') as f:
    f.write(f90code(-A, 'A') + "\n")
    f.write(f90code(A.jacobian([rho,T,H]), 'M') + "\n")
