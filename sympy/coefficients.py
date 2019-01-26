#!/usr/bin/env python
# coding: utf-8

from sympy import Symbol, Function, Lambda, Derivative, IndexedBase, Eq, symbols, Integer, Rational, Matrix, MatrixSymbol, Wild, simplify, sqrt, exp, log, Piecewise
from numpy import ndarray, zeros, linspace, logspace, meshgrid
from sys import stdin, stdout, stderr
from io import StringIO

#------------------------------------------------------------------------------#

# Real Positive Symbol
RPS = lambda x: Symbol(x, real = True, positive = True)

alpha = RPS('alpha')
eta   = RPS('eta')
nu    = RPS('nu')
omega = RPS('omega')
F_acc = RPS('F_acc')
mu    = RPS('mu')
m_el  = RPS('m_el')
m_H   = RPS('m_H')
c     = RPS('c')
sigma = RPS('sigma')
radius = RPS('radius')
rschw = RPS('rschw')
kboltz     = RPS('k_B')
use_qmri = Symbol('use_quench_mri')
use_prad_in_alpha = Symbol('use_prad_in_alpha')
use_relcompt = Symbol('use_precise_balance')
use_fluxcorr = Symbol('use_flux_correction')

global_variables = {
    F_acc:  'facc',
    m_el:   'cgs_mel',
    m_H:    'cgs_mhydr',
    c:      'cgs_c',
    sigma:  'cgs_stef',
    kboltz: 'cgs_boltz',
    mu:     'miu',
    alpha:  'alpha',
    eta:    'zeta',
}

global_expressions = [
    (kboltz / ( m_el * c**2 ), RPS('cgs_k_over_mec2')),
    (kboltz / ( m_el * c ), RPS('cgs_k_over_mec')),
    (kboltz / m_H, RPS('cgs_k_over_mh')),
]

#------------------------------------------------------------------------------#

from sympy.printing.fcode import FCodePrinter
from sympy.printing.precedence import precedence, precedence

# this class is to fix some issues with FCodePrinter from sympy
class F90CodePrinter(FCodePrinter):
    def __init__(self, settings = {}):
        settings['source_format'] = 'free'
        settings['standard'] = 2008
        super(F90CodePrinter, self).__init__(settings)

    def _print_MatrixElement(self, expr):
        # if we have a vector (second dimension = 1), print only one index
        if expr.parent.shape[1] <= 1:
            return "{0}({1})".format(expr.parent, expr.i + 1)
        return "{0}({1},{2})".format(expr.parent, expr.i + 1, expr.j + 1)

    def _print_Piecewise(self, expr):
        from sympy import RealNumber
        # this fix is to prevent sympy from printing constructs like
        # merge(x,0,condition) which are illegal in Fortran
        newargs = ((RealNumber(x) if x.is_Integer else x \
            for x in y) for y in expr.args)
        return super(F90CodePrinter, self)._print_Piecewise(Piecewise(*newargs))

fprinter = F90CodePrinter()

#------------------------------------------------------------------------------#

fcoeff = open('../src/mrxcoeff.fi','w')

fsub_coeff = """pure subroutine mrx_coeff1_{name} (z,Y,D,F,A,MY,MD)
use iso_fortran_env, only: dp => real64
implicit none
real(dp), intent(in) :: z
! first column: values, second column: 1ord derivatives
real(dp), dimension(:), intent(in) :: Y,D
! each row is one Function: its value and derivatives to rho and T
real(dp), dimension(:,:), intent(in) :: F
! right-hand-side of the equation
real(dp), dimension(:), intent(out) :: A
! jacobians with respect to Y and dY/dz
real(dp), dimension(:,:), intent(out) :: MY, MD
{A}\n{MY}\n{MD}\nend subroutine\n\n"""

fsub_bound = """pure subroutine mrx_coeff{lr}_{name} (z,Y,F,B,MB)
use iso_fortran_env, only: dp => real64
implicit none
real(dp), intent(in) :: z
real(dp), dimension(:), intent(in) :: Y
real(dp), dimension(:,:), intent(in) :: F
real(dp), dimension(:), intent(out) :: B
real(dp), dimension(:,:), intent(out) :: MB
{B}\n{MB}\nend subroutine\n\n"""

fsub_constr = """pure subroutine mrx_coeff0_{name} (z,Y,F,C,MC)
use iso_fortran_env, only: dp => real64
implicit none
real(dp), intent(in) :: z
real(dp), dimension(:), intent(in) :: Y
real(dp), dimension(:,:), intent(in) :: F
real(dp), dimension(:), intent(out) :: C
real(dp), dimension(:,:), intent(out) :: MC
{B}\n{MB}\nend subroutine\n\n"""

fsub_yout = """pure subroutine mrx_output_{name} (z,Y,YY)
use iso_fortran_env, only: dp => real64
implicit none
real(dp), intent(in) :: z
real(dp), dimension(:), intent(in) :: Y
real(dp), dimension(:), intent(out) :: YY
{YY}\nend subroutine\n\n"""

#------------------------------------------------------------------------------#

fswptrs = open('../src/mrxptrs.fi','w')
fswdims = open('../src/mrxdims.fi','w')
fswhash = open('../src/mrxhash.fi','w')
fswfout = open('../src/mrxfout.fi','w')
fswall = [fswptrs, fswdims, fswhash, fswfout]

for f in fswall: f.write("select case (nr)\n")

#------------------------------------------------------------------------------#

choices = [
#   BIL     FULL   MAGN   CND
    (False, False, False, False),
    (False, False, True,  False),
    (True,  False, True,  False),
    (True,  True,  True,  False),
    # (False, False, False, True ),
    # (True,  True,  True,  True ),
    # (True,  False, True,  True ),
]

for balance, bilfull, magnetic, conduction in choices:

    if bilfull and not balance: continue

    # współrzędne, po której liczymy
    z       = Symbol('z')

    # niewiadome
    rho     = Function('rho')(z)
    T_gas   = Function('T')(z)
    F_rad   = Function('F_rad')(z)

    T_rad   = Function('T_rad')(z) if balance else T_gas
    P_mag   = Function('P_mag')(z) if magnetic else 0
    F_cond  = Function('F_cond')(z) if conduction else 0

    #--------------------------------------------------------------------------#

    # niewiadome w postaci listy
    yvar = [ rho, T_gas ]
    if balance: yvar.append(T_rad)
    yvar.append(F_rad)
    if conduction: yvar.append(F_cond)
    if magnetic: yvar.append(P_mag)

    yvar_all = [ rho, T_gas, T_rad, F_rad, P_mag, F_cond ]
    yvar_lbl = ['rho', 'temp', 'trad', 'frad', 'pmag', 'fcnd']
    yval_hash = [ yvar.index(y)+1 if y in yvar else 0 for y in yvar_all ]

    #--------------------------------------------------------------------------#

    # nieprzezroczystosci i rpzewodnizctwa
    kabs = Function('fkabs')(rho,T_gas)
    kabp = Function('fkabp')(rho,T_gas)
    ksct = Function('fksct')(rho,T_gas)
    kcnd = Function('fkcnd')(rho,T_gas)

    fvec = [ kabs, kabp, ksct, kcnd ]

    fval = MatrixSymbol('F', 3, len(fvec))

    def Function_derivatives(eq):
        Functions_extra = []
        for ifun,f in zip(range(len(fvec)),fvec):
            neq1rgs = len(f.args)
            ff = f.func
            w = [ Wild('w{}'.format(i)) for i in range(neq1rgs) ]
            for i in range(neq1rgs):
                replace_what = Derivative(ff(*w),w[i])
                eq = eq.replace( replace_what, fval[i+1,ifun] )
            eq = eq.replace( ff(*w), fval[0,ifun] )
        return eq.doit()

    #--------------------------------------------------------------------------#

    # predkosc wznoszenia pola
    vrise = Lambda(z, eta * omega * z)
    # cisnienie gazu
    P_gas = kboltz * rho * T_gas / ( mu * m_H )
    # predkosc dzwieko
    csound = sqrt(P_gas / rho)
    # cisnienie promieniowania
    P_rad = 4 * sigma / (3 * c) * T_rad**4
    # plasma beta
    beta = P_gas / P_mag
    # radiative beta
    betarad = P_gas / P_rad
    # strumien Poyntinga
    F_mag = 2 * P_mag * vrise(z)
    # strumien calkowity
    F_tot = F_rad + F_mag + F_cond

    # cisnienie calkowite w alpha-prescription
    P_tot_gen = P_gas + P_mag \
        + Piecewise((P_rad, use_prad_in_alpha), (0.0, True))

    betamri = 2 * csound / (omega * radius * rschw)
    qmri0 = 1 / (1 + (betamri / beta)**4)
    qmri = Piecewise((qmri0, use_qmri), (1.0, True))

    #--------------------------------------------------------------------------#

    eq1ord = []
    eq0ord = []
    boundL = []
    boundR = []

    #--------------------------------------------------------------------------#
    # Hydrostatic equilibrium
    #
    eq1ord.append(P_gas.diff(z)      \
            - kabs * rho / c * F_rad   \
            + Derivative(P_mag,z)       \
            + omega**2 * rho * z)

    #--------------------------------------------------------------------------#
    # Dyfuzja promieniowania
    #
    eq1ord.append(
        3 * kabs * rho * F_rad + 16 * sigma * T_rad**3 * Derivative(T_rad,z)
    )

    #--------------------------------------------------------------------------#
    # Bilans grzania i chlodzenia
    #
    if magnetic:
        heat = 2 * P_mag * vrise(z).diff(z) \
            - alpha * omega * (P_tot_gen * qmri - 2 * nu * P_mag)
    else:
        heat = alpha * omega * P_tot_gen

    eq1ord.append(
        Derivative(F_rad,z) + Derivative(F_cond,z) - heat
    )

    q = 2 + (alpha / eta) * (nu - 1)
    # if high in the atmosphere H(z) ~ 1 / z**q, then missing flux from
    # ztop to infinity yields Fcorr = H(ztop) * ztop / (q - 1)
    F_tot_fix = Piecewise(((alpha * omega * P_mag) * z / (q - 1), use_fluxcorr), (0.0, True))
    F_rad_fix = Piecewise((heat * z / (q - 1), use_fluxcorr), (0.0, True))

    boundL.append(F_rad)
    boundR.append(F_tot + (F_tot_fix if magnetic else 0)  - F_acc)
    boundR.append(F_rad - 2 * sigma * T_rad**4)


    #--------------------------------------------------------------------------#
    # Magnetic field evolution
    #
    if magnetic:
        eq1ord.append(
            2 * P_mag * vrise(z).diff(z)        \
            + vrise(z) * Derivative(P_mag,z)    \
            - alpha * omega * (P_tot_gen * qmri - nu * P_mag)
        )
        boundL.append(  2 * P_mag * vrise(z).diff(z)            \
                        - alpha * omega * (P_tot_gen - nu * P_mag)  )

    #--------------------------------------------------------------------------#
    # Funkcja zrodlowa
    #
    if balance:
        sdyf = kabp * (T_gas**4 - T_rad**4)
        relcor = 1 + 4 * kboltz * T_gas / (m_el * c**2)
        relcor = Piecewise((relcor, use_relcompt), (1.0, True))
        ssct = ksct * T_rad**4 * kboltz / (m_el * c**2) \
            * 4 * (relcor * T_gas - T_rad)
        radcool = 4 * sigma * rho * ((sdyf if bilfull else 0) + ssct)

        if conduction:
            eq1ord.append(radcool - Derivative(F_rad,z))
            boundL.append(T_gas - T_rad)
        else:
            eq0ord.append(heat - radcool)

    #--------------------------------------------------------------------------#
    # Dyfuzja ciepla przez przewodnictwo
    #
    if conduction:
        eq1ord.append(
            F_cond + kcnd * Derivative(T_gas,z)
        )
        boundL.append(F_cond)

    #--------------------------------------------------------------------------#

    neq1 = len(eq1ord)
    neq0 = len(eq0ord)
    ny = len(yvar)
    nbl = len(boundL)
    nbr = len(boundR)

    assert neq1 + neq0 == ny
    assert neq1 == nbl + nbr

    Y = MatrixSymbol('Y', ny, 1)
    D = MatrixSymbol('D', ny, 1)

    #--------------------------------------------------------------------------#

    def discretize(eq):
        if eq == 0: return eq
        for iy,y in zip(range(ny),yvar):
            eq = eq.subs(Derivative(y,z),D[iy]).subs(y,Y[iy])
        for k,v in global_variables.items():
            eq = eq.subs(k,Symbol(v, real = True, positive = k.is_positive))
        for k,v in global_expressions:
            eq = eq.subs(k,v)
        return Function_derivatives(eq)

    #--------------------------------------------------------------------------#

    model_name = '{magn}{comp}{cond}'.format(
        magn = 'm' if magnetic else 'a',
        comp = ('c' if bilfull else 'w') if balance else 'd',
        cond = 't' if conduction else '',
    )

    model_nr = 1 \
        + (1 if balance else 0) \
        + (2 if bilfull else 0) \
        + (4 if magnetic else 0) \
        + (8 if conduction else 0)
    print("{:4d} -> {}".format(model_nr,model_name.upper()))

    #--------------------------------------------------------------------------#

    fcoeff.write('!' + 78 * '-' + '!\n')
    fcoeff.write("! {}\n".format("heating-cooling balance ({})".format("full equation" if bilfull else "compton term only") if balance else "thermal diffusion"))
    fcoeff.write("! {}\n".format("magnetic heating" if magnetic else "alpha prescription"))
    fcoeff.write("! {}\n".format("radiation + thermal conduction" if conduction \
        else "no thermal conduction"))

    #--------------------------------------------------------------------------#

    A = Matrix(zeros((neq1,)))
    MY = Matrix(zeros((neq1,ny)))
    MD = Matrix(zeros((neq1,ny)))

    for ieq,eq in zip(range(neq1),eq1ord):
        A[ieq] = discretize(-eq)
        for iy,y in zip(range(ny),yvar):
            MY[ieq,iy] = discretize(eq.diff(y))
            MD[ieq,iy] = discretize(eq.diff(Derivative(y,z)))

    fcoeff.write(fsub_coeff.format(
        name = model_name,
        A = fprinter.doprint(A, 'A'),
        MY = fprinter.doprint(MY, 'MY'),
        MD = fprinter.doprint(MD, 'MD'),
    ))

    #--------------------------------------------------------------------------#

    if neq0 > 0:
        C = Matrix(zeros((neq0,)))
        MC = Matrix(zeros((neq0,ny)))

        for ict,ct in zip(range(neq0),eq0ord):
            C[ict] = discretize(-ct)
            for iy,y in zip(range(ny),yvar):
                MC[ict,iy] = discretize(ct.diff(y))

        fcoeff.write(fsub_constr.format(
            name = model_name,
            B = fprinter.doprint(C, 'C'),
            MB = fprinter.doprint(MC, 'MC'),
        ))

    #--------------------------------------------------------------------------#

    if nbl > 0:
        BL = Matrix(zeros((nbl,)))
        MBL = Matrix(zeros((nbl,ny)))

        for ib,b in zip(range(nbl),boundL):
            BL[ib] = discretize(-b)
            for iy,y in zip(range(ny),yvar):
                MBL[ib,iy] = discretize(b.diff(y))

        fcoeff.write(fsub_bound.format(
            name = model_name, lr = 'bl',
            B = fprinter.doprint(BL, 'B'),
            MB = fprinter.doprint(MBL, 'MB'),
        ))

    #--------------------------------------------------------------------------#

    if nbr > 0:
        BR = Matrix(zeros((nbr,)))
        MBR = Matrix(zeros((nbr,ny)))

        for ib,b in zip(range(nbr),boundR):
            BR[ib] = discretize(-b)
            for iy,y in zip(range(ny),yvar):
                MBR[ib,iy] = discretize(b.diff(y))

        fcoeff.write(fsub_bound.format(
            name = model_name, lr = 'br',
            B = fprinter.doprint(BR, 'B'),
            MB = fprinter.doprint(MBR, 'MB'),
        ))

    #--------------------------------------------------------------------------#

    yout = [
        rho,  T_gas, T_rad,
        P_gas, P_rad, P_mag,
        F_rad, F_mag, F_cond,
        P_tot_gen, heat, vrise(z), qmri0 if magnetic else 0,
    ]

    fcoeff.write(fsub_yout.format(
        name = model_name,
        YY = fprinter.doprint(Matrix([ discretize(y) for y in yout ]), 'YY'),
    ))

    #--------------------------------------------------------------------------#

    for f in fswall: f.write('case({})\n'.format(model_nr))

    fswptrs.write ('  feq0 => {}\n'      \
        .format('mrx_coeff0_' + model_name if neq0 > 0 else 'NULL()'))
    fswptrs.write ('  feq1 => {}\n'.format('mrx_coeff1_' + model_name))
    fswptrs.write ('  fbl  => {}\n'      \
        .format('mrx_coeffbl_' + model_name if nbl > 0 else 'NULL()'))
    fswptrs.write ('  fbr  => {}\n'      \
        .format('mrx_coeffbr_' + model_name if nbr > 0 else 'NULL()'))
    fswfout.write('  fout => {}\n'.format('mrx_output_' + model_name))

    fswdims.write ('  ny   = {}\n'.format(ny))
    fswdims.write ('  neq0 = {}\n'.format(neq0))
    fswdims.write ('  neq1 = {}\n'.format(neq1))
    fswdims.write ('  nbl  = {}\n'.format(nbl))
    fswdims.write ('  nbr  = {}\n'.format(nbr))

    fswhash.write ('  ihash(:) = [{}]\n'.format(", ".join([ str(i) \
            for i in yval_hash ])))

#------------------------------------------------------------------------------#

for f in fswall:
    f.write("case default\n  error stop \"this model is not available\"\nend select\n")
    f.close()
