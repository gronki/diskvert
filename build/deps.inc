alphadisk.o: balance.o cgs.o fileunits.o globals.o grid.o odeintegr.o rk4settings.o ss73solution.o thrfun.o
balance.o: cgs.o globals.o
balance.o: globals.o thrfun.o
diskvert.o: balance.o fileunits.o globals.o grid.o integrderiv.o interpol.o relaxation.o settings.o ss73solution.o
dv-alpha.o: alphadisk.o fileunits.o globals.o integrderiv.o results.o settings.o
dv-mag.o: fileunits.o globals.o grid.o modelmag.o results.o settings.o
dv-rad1.o: alphadisk.o cgs.o fileunits.o globals.o settings.o ss73solution.o
globals.o: cgs.o
modelmag.o: balance.o cgs.o globals.o odeintegr.o rk4settings.o settings.o ss73solution.o thrfun.o
relaxation.o: cgs.o globals.o mrxcoeff.fi mrxdims.fi mrxfout.fi mrxhash.fi mrxptrs.fi thrfun.o
rk4settings.o: globals.o
settings.o: fileunits.o globals.o
ss73solution.o: cgs.o globals.o

diskvert: balance.o cgs.o diskvert.o fileunits.o globals.o grid.o integrderiv.o interpol.o relaxation.o settings.o ss73solution.o thrfun.o
dv-alpha: alphadisk.o balance.o cgs.o dv-alpha.o fileunits.o globals.o grid.o integrderiv.o odeintegr.o results.o rk4settings.o settings.o ss73solution.o thrfun.o
dv-mag: balance.o cgs.o dv-mag.o fileunits.o globals.o grid.o modelmag.o odeintegr.o results.o rk4settings.o settings.o ss73solution.o thrfun.o
dv-rad1: alphadisk.o balance.o cgs.o dv-rad1.o fileunits.o globals.o grid.o odeintegr.o rk4settings.o settings.o ss73solution.o thrfun.o
