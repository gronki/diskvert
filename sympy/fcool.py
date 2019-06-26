from sympy import *

var('cgs_c cgs_kapes cgs_k_over_mec2 cgs_stef', real = True, positive = True)
var('rho temp trad heat', real = True, positive = True)
var('use_relcompt')

kabp = Function('kabp')(rho,temp)
ksct = Function('ksct')(rho,temp)

yyb = kabp * (temp**4 - trad**4)
relcor = Piecewise((1 + 4 * cgs_k_over_mec2 * temp, use_relcompt), (1.0, True))
yyc = ksct * trad**4 * cgs_k_over_mec2 * 4 * (temp * relcor - trad)
yy = 4 * cgs_stef * rho * (yyb + yyc)

kabpv = IndexedBase('kabpv')
ksctv = IndexedBase('ksctv')

def ff(x,y):
    x = x.subs(Derivative(ksct, rho),  ksctv[2])
    x = x.subs(Derivative(ksct, temp), ksctv[3])
    x = x.subs(ksct, ksctv[1])
    x = x.subs(Derivative(kabp, rho),  kabpv[2])
    x = x.subs(Derivative(kabp, temp), kabpv[3])
    x = x.subs(kabp, kabpv[1])
    return fcode(simplify(x), assign_to = y, source_format = 'free',
            standard = 2008, contract = False)

with open('fcool.f90', 'w') as f:
    f.write(ff(yy,            'cool'      ) + '\n')
    f.write(ff(yy.diff(rho),  'cool_drho' ) + '\n')
    f.write(ff(yy.diff(temp), 'cool_dtemp') + '\n')
