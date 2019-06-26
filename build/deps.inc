modelmag.o: balance.o ss73solution.o rk4settings.o thrfun.o odeintegr.o settings.o cgs.o globals.o
dv-alpha.o: alphadisk.o fileunits.o results.o settings.o globals.o
alphadisk.o: balance.o ss73solution.o rk4settings.o thrfun.o odeintegr.o fileunits.o cgs.o grid.o globals.o
diskvert.o: balance.o integrderiv.o ss73solution.o fileunits.o settings.o grid.o interpol.o globals.o relaxation.o
balance.o: globals.o thrfun.o
globals.o: cgs.o
dv-mag.o: fileunits.o results.o settings.o grid.o globals.o modelmag.o
settings.o: fileunits.o globals.o
rk4settings.o: globals.o
ss73solution.o: cgs.o globals.o
dv-rad1.o: ss73solution.o alphadisk.o fileunits.o settings.o cgs.o globals.o
balance.o: cgs.o globals.o
relaxation.o: mrxhash.fi thrfun.o mrxcoeff.fi mrxfout.fi cgs.o mrxptrs.fi globals.o mrxdims.fi
