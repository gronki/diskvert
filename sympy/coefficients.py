#!/usr/bin/env python
# coding: utf-8

from sympy import Symbol, Function, Lambda, Derivative, IndexedBase, Eq, symbols, \
    Integer, Rational, Matrix, MatrixSymbol, Wild, simplify, sqrt, exp, log, Piecewise, var
from numpy import ndarray, zeros, linspace, logspace, meshgrid
from sys import stdin, stdout, stderr
from io import StringIO

# ------------------------------------------------------------------------------#

# these variable names correspond (mostly) to names in Fortran code
var('alpha eta nu mbh mdot radius rschw omega facc miu', real=True, positive=True)
var('cgs_mel cgs_mhydr cgs_c cgs_stef cgs_boltz', real=True, positive=True)
var('qmri_kill zbreak threshpow condux hydrox', real=True, positive=True)
var('cgs_k_over_mec2 cgs_k_over_mec cgs_k_over_mh', real=True, positive=True)
var('use_prad_in_alpha use_relcompt use_quench_mri use_flux_correction')

# ------------------------------------------------------------------------------#

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

# ------------------------------------------------------------------------------#

# headers of Fortran procedures
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

# ------------------------------------------------------------------------------#

# .fi files are outputs of this script

# .fi files contain the code generated with this script
# you have to regenerate them only if you want to change
# anyhting in the equations

# mrxcoeff.fi will contain the partial derivatives of all equations
# the option suffixes are explained in the output file
fcoeff = open('../src/mrxcoeff.fi','w')
# this will choose procedures to compute matrix coefficients
# for boundary conditions and points inside the interval
fswptrs = open('../src/mrxptrs.fi','w')
# this file will select the number of variables (ny),
# differental equations of N-th order (neqN)
# equat. plane (nbl) and surface (nbr) boundary conditions
fswdims = open('../src/mrxdims.fi','w')
# ...
fswhash = open('../src/mrxhash.fi','w')
# this file makes the choice of procedure to be used for
fswfout = open('../src/mrxfout.fi','w') 
fswall = [fswptrs, fswdims, fswhash, fswfout]

for f in fswall: f.write("select case (nr)\n")

# ------------------------------------------------------------------------------#

# BIL - if we have two temperatures (Trad, T) or only one
# FULL - if we use proper equation for Lambda or simplified one
# MAGN - if we have magnetic disk or only alpha prescription
# CND - if we include thermal conduction - does not work
# MASSFLX - mass flux
# the list below is a list of model possibilities we generate
choices = [
    #   BIL     FULL   MAGN   CND    MASSFLX
    (False, False, False, False, False),
    (False, False, True, False,  False),  # default
    (False, False, True, False,  True),
    (True,  False, True, False,  False),
    (True,  False, True, False,  True),
    (True,  True,  True, False,  False),
    (True,  True,  True, False,  True),
    (False, False, True, True,   False),
    (True,  True,  True, True,   False),
]

for balance, bilfull, magnetic, conduction, massflux in choices:

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

    # --------------------------------------------------------------------------#

    # niewiadome w postaci listy
    yvar = [ rho, T_gas ]
    if balance: yvar.append(T_rad)
    yvar.append(F_rad)
    if conduction: yvar.append(F_cond)
    if magnetic: yvar.append(P_mag)

    yvar_all = [ rho, T_gas, T_rad, F_rad, P_mag, F_cond ]
    yvar_lbl = ['rho', 'temp', 'trad', 'frad', 'pmag', 'fcnd']
    yval_hash = [ yvar.index(y)+1 if y in yvar else 0 for y in yvar_all ]

    # --------------------------------------------------------------------------#

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

    # --------------------------------------------------------------------------#

    # cisnienie gazu
    P_gas = (cgs_k_over_mh / miu) * (rho * T_gas)
    # cisnienie promieniowania
    P_rad = (4 * cgs_stef) / (3 * cgs_c) * T_rad**4
    # plasma beta
    beta = P_gas / P_mag
    # radiative beta
    betarad = P_gas / P_rad
    # predkosc wyplywu
    vmag = lambda z: eta * omega * z
    # vmag = lambda z: eta * omega * Piecewise((z - zbreak, z >= zbreak), (0., True))
    # strumien Poyntinga
    F_mag = 2 * P_mag * vmag(z)
    # strumien calkowity
    F_tot = F_rad + F_mag + F_cond

    # cisnienie calkowite w alpha-prescription
    P_therm = P_gas + Piecewise((P_rad, use_prad_in_alpha), (0, True))
    P_gen = P_therm + (P_mag if magnetic else 0)

    # if x = hi/lo > 1 then thresh = 1 else thresh = 0
    # thresh = lambda x: x**threshpow / (1 + x**threshpow)
    thresh = lambda x: (x + x**threshpow) / (2 + x + x**threshpow)

    # betamri = 2 * cs / vr
    P_mag_max = sqrt(1.67 * P_gas * rho) * (omega * radius * rschw) / 2
    qmri0 = thresh(P_mag_max / P_mag)
    qmri = Piecewise((qmri0 * qmri_kill + (1.0 - qmri_kill), use_quench_mri), (1.0, True))

    # variable parameters
    valpha = qmri * alpha
    veta = vmag(z).diff(z) / omega
    vnu = nu

    # --------------------------------------------------------------------------#

    eq1ord = []
    eq0ord = []
    boundL = []
    boundR = []

    # --------------------------------------------------------------------------#
    # Hydrostatic equilibrium
    #
    # dtrad_dz = Derivative(T_rad, z)
    dtrad_dz = - 3 * kabs * rho * F_rad / (16 * cgs_stef * T_rad**3) 
    # dT_dz = Piecewise((dtrad_dz, use_simple_hydro), (Derivative(T_gas, z), True))
    dT_dz = dtrad_dz * (1 - hydrox) + Derivative(T_gas, z) * hydrox
    pgas_gradient = (cgs_k_over_mh / miu) * (Derivative(rho, z) * T_gas + rho * dT_dz)
    eq1ord.append(pgas_gradient
            - kabs * rho / cgs_c * F_rad   \
            + Derivative(P_mag, z)     \
            + omega**2 * rho * z)

    # --------------------------------------------------------------------------#
    # Dyfuzja promieniowania
    #
    eq1ord.append(
        3 * kabs * rho * F_rad + 16 * cgs_stef * T_rad**3 * Derivative(T_rad,z)
    )

    # --------------------------------------------------------------------------#
    # Magnetic field evolution
    #
    if magnetic:
        etazdpmag = (2 * veta + alpha * vnu) * P_mag - valpha * P_gen
        eq1ord.append(etazdpmag + (vmag(z) / omega) * Derivative(P_mag, z))
        boundL.append(etazdpmag)

    # --------------------------------------------------------------------------#
    # Bilans grzania i chlodzenia
    #
    if magnetic:
        heatm = (2 * veta + alpha * vnu) * omega * P_mag - valpha * omega * P_gen
        heatr = alpha * vnu * omega * P_mag 
        heat = heatm + heatr
    else:
        heat = alpha * omega * P_therm

    if not (conduction and balance): eq1ord.append(Derivative(F_rad,z) + Derivative(F_cond,z) - heat)

    q  = 2 + (alpha / eta) * (nu - 1)
    q1 = 2 + (alpha / eta) * nu
    # q1 = Piecewise((q1, use_quench_mri), (q, True))
    # if high in the atmosphere H(z) ~ 1 / z**q, then missing flux from
    # ztop to infinity yields Fcorr = H(ztop) * ztop / (q - 1)
    F_tot_fix = Piecewise(((alpha * omega * P_mag) * z / (q1 - 1), use_flux_correction), (0.0, True))
    F_rad_fix = Piecewise((heat * z / (q1 - 1), use_flux_correction), (0.0, True))

    boundL.append(F_rad)
    boundR.append(F_tot + (F_tot_fix if magnetic else 0)  - facc)
    boundR.append(F_rad - 2 * cgs_stef * T_rad**4)

    # --------------------------------------------------------------------------#
    # Funkcja zrodlowa
    #
    if balance:
        sdyf = kabp * (T_gas**4 - T_rad**4)
        relcor = 1 + 4 * cgs_k_over_mec2 * T_gas
        relcor = Piecewise((relcor, use_relcompt), (1.0, True))
        ssct = ksct * T_rad**4 * 4 * cgs_k_over_mec2 * (relcor * T_gas - T_rad)
        radcool = 4 * cgs_stef * rho * ((sdyf if bilfull else 0) + ssct)

        if conduction:
            eq1ord.append(radcool + Derivative(F_cond,z) - heat) #############
            eq1ord.append(radcool - Derivative(F_rad,z))
            boundL.append(heat - radcool)
        else:
            eq0ord.append(heat - radcool)

    # --------------------------------------------------------------------------#
    # Dyfuzja ciepla przez przewodnictwo
    #
    if conduction:
        eq1ord.append(
            F_cond + kcnd * (Derivative(T_rad,z) * (1 - condux) + Derivative(T_gas,z) * condux)
        )
        boundL.append(F_cond)

    # --------------------------------------------------------------------------#

    neq1 = len(eq1ord)
    neq0 = len(eq0ord)
    ny = len(yvar)
    nbl = len(boundL)
    nbr = len(boundR)

    assert neq1 + neq0 == ny
    assert neq1 == nbl + nbr

    Y = MatrixSymbol('Y', ny, 1)
    D = MatrixSymbol('D', ny, 1)

    # --------------------------------------------------------------------------#

    def discretize(eq):
        if eq == 0: return eq
        for iy,y in zip(range(ny),yvar):
            eq = eq.subs(Derivative(y,z),D[iy]).subs(y,Y[iy])
        return Function_derivatives(eq)

    # --------------------------------------------------------------------------#

    model_name = "{magn}{comp}{cond}{flx}".format(
        magn="m" if magnetic else "a",
        comp=("c" if bilfull else "w") if balance else "d",
        cond="t" if conduction else "",
        flx="f" if massflux else "",
    )

    model_nr = (
        1
        + (1 if balance else 0)
        + (2 if bilfull else 0)
        + (4 if magnetic else 0)
        + (8 if conduction else 0)
        + (16 if massflux else 0)
    )
    print("{:4d} -> {}".format(model_nr,model_name.upper()))

    # --------------------------------------------------------------------------#

    fcoeff.write('!' + 78 * '-' + '!\n')
    fcoeff.write("! {}\n".format("heating-cooling balance ({})".format("full equation" if bilfull else "compton term only") if balance else "thermal diffusion"))
    fcoeff.write("! {}\n".format("magnetic heating" if magnetic else "alpha prescription"))
    fcoeff.write("! {}\n".format("radiation + thermal conduction" if conduction \
        else "no thermal conduction"))

    # --------------------------------------------------------------------------#

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

    # --------------------------------------------------------------------------#

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

    # --------------------------------------------------------------------------#

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

    # --------------------------------------------------------------------------#

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

    # --------------------------------------------------------------------------#

    yout = [
        rho,  T_gas, T_rad,
        P_gas, P_rad, P_mag,
        F_rad, F_mag, F_cond,
        P_gen, heat, heatm if magnetic else 0, heatr if magnetic else 0, 
        vmag(z) if magnetic else 0, 
        qmri0 if magnetic else 0, 
    ]

    fcoeff.write(fsub_yout.format(
        name = model_name,
        YY = fprinter.doprint(Matrix([ discretize(y) for y in yout ]), 'YY'),
    ))

    # --------------------------------------------------------------------------#

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

# ------------------------------------------------------------------------------#

for f in fswall:
    f.write("case default\n  error stop \"this model is not available\"\nend select\n")
    f.close()
