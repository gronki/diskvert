import numpy as np
import numpy.linalg
from sympy import symbols, Symbol, Rational, Number, log as sympy_log, \
    expand_log, Matrix, exp as sympy_exp, nsimplify, solveset, fcode, Eq, \
    sqrt as sympy_sqrt, Function, var
from diskvert.cgs import *

mbh, mdot, r, alpha = symbols('mbh mdot r alpha', real = True, positive = True)

logY_all = list()
eqstream = list()

var('cgs_kapes kappa_abs_0', real = True, positive = True)
kappa_es = cgs_kapes

ff = Function('f')

for region in [1,2,3]:
    rho, T, H = symbols('rho{n} T{n} H{n}'.format(n = region), real = True, positive = True)

    r_cgs = r * mbh * sol_rschw
    mdot_cgs = mdot * mbh * sol_mdot_edd
    omega = ( cgs_graw * cgs_msun * mbh / r_cgs**3 )**0.5
    f = Symbol('f', real = True, positive = True)

    yvar = [rho, T, H]
    log_yvar = [ sympy_log(x) for x in yvar ]

    P1 = 4 * cgs_stef / (3 * cgs_c) * T**4
    P23 = 2 * cgs_boltz * rho * T / cgs_mhydr
    kap12 = kappa_es
    kap3 = kappa_abs_0 * rho * T**Rational(-7,2)

    # in regions 1,2 electron scattering opacity, else free-free
    kap = kap3 if region == 3 else kap12
    # radiation pressure in the inner region, gas press everywhere else
    P =     P1 if region == 1 else P23

    # left and right hand sides of 3 equations
    L1 = 4 * cgs_pi * alpha * rho * H**3 * omega
    R1 = mdot_cgs * f
    L2 = ( Rational(3,4) * rho * H**2 )**2 * alpha * omega **3 * kap
    R2 = 4 * cgs_stef * T**4
    L3 = P
    R3 = rho * H**2 * omega**2

    # obtain so that RHS = 1 everywhere and we can logarithm
    A = [ L1/R1, L2/R2, L3/R3 ]
    # logarithm the lhs
    logA = [ expand_log(sympy_log(x)) for x in A ]

    # now collect the coefficients next to logarithms of each
    # variable: rho, T, H
    M = np.ndarray((3,3))
    for i in range(3):
        for j in range(3):
            M[i,j] = logA[i].coeff(log_yvar[j])
    print M

    # the above will constitute our new left hand sides
    lhs = [ sum([ M[i,j] * log_yvar[j] for j in range(3) ]) for i in range(3) ]
    # right hand sides shall be a remainder
    rhs = [ lhs[i] - logA[i] for i in range(3) ]

    for i in range(3):
        print "{} = {}".format(lhs[i],rhs[i])

    # invert the matrix
    Minv = numpy.linalg.inv(M)
    # simplify the inversion coefficients before multiplication, e.g. 0.2 -> 1/5
    # calculate the variables
    logY = Matrix([ [ nsimplify(Minv[i,j]) for j in range(3) ] for i in range(3) ]) * Matrix(rhs)
    # add to the global array. we will need them later
    logY_all.append([logY[i] for i in range(3)])

    for i in range(3):
        eqstream.append(Eq(yvar[i], sympy_exp(logY[i].subs(f,ff(r)))))

temp = Symbol('logR')

r12,r23 = symbols('r12 r23', real = True, positive = True)
f12,f23 = symbols('f12 f23', real = True, negative = False)

# przejscie 1->2
eq = nsimplify(logY_all[1][0]) - 3 * nsimplify(logY_all[1][1]) - sympy_log( (2 * cgs_stef * cgs_mhydr) / (3 * cgs_boltz * cgs_c) )
for x in solveset(eq.replace(sympy_log(r), temp),temp):
    eqstream.append(Eq(r12, sympy_exp(x.subs(f,ff(r12)))))

# przejscie 2->3
eq = nsimplify(logY_all[1][0]) - Rational(7,2) * nsimplify(logY_all[1][1]) - sympy_log( kappa_es / kappa_abs_0 )
for x in solveset(eq.replace(sympy_log(r), temp),temp):
    eqstream.append(Eq(r23, sympy_exp(x.subs(f,ff(r23)))))

with open('3zones.f90','w') as f:
    for eq in eqstream:
        f.write(fcode(eq.rhs.evalf(6), eq.lhs, standard = 2008,
            source_format = 'free', contract = False, user_functions = {'f': 'f'}) + "\n")
