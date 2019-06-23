ss73solution.o: globals.o cgs.o
diskvert.o: balance.o settings.o interpol.o grid.o fileunits.o ss73solution.o globals.o integrderiv.o relaxation.o
dv-mag.o: settings.o grid.o fileunits.o modelmag.o globals.o results.o
balance.o: globals.o cgs.o
globals.o: cgs.o
dv-alpha.o: settings.o fileunits.o alphadisk.o globals.o results.o
balance.o: globals.o thrfun.o
modelmag.o: cgs.o balance.o settings.o ss73solution.o rk4settings.o odeintegr.o globals.o thrfun.o
dv-alpha-rx.o: settings.o grid.o fileunits.o ss73solution.o globals.o relaxation.o
settings.o: globals.o fileunits.o
dv-rad1.o: cgs.o settings.o fileunits.o ss73solution.o alphadisk.o globals.o
alphadisk.o: cgs.o balance.o grid.o fileunits.o ss73solution.o rk4settings.o odeintegr.o globals.o thrfun.o
rk4settings.o: globals.o
relaxation.o: mrxdims.fi cgs.o mrxhash.fi mrxcoeff.fi mrxptrs.fi globals.o mrxfout.fi thrfun.o
