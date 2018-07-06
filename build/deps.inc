dv-mag-rx.o: settings.o globals.o balance.o ss73solution.o fileunits.o         \
    relaxation.o grid.o integrderiv.o interpol.o 
balance.o: globals.o thrfun.o 
dv-mag.o: settings.o globals.o modelmag.o fileunits.o grid.o results.o 
dv-rad1.o: settings.o globals.o cgs.o alphadisk.o ss73solution.o fileunits.o 
balance.o: globals.o cgs.o 
ss73solution.o: globals.o cgs.o 
dv-alpha-rx.o: settings.o globals.o ss73solution.o fileunits.o relaxation.o    \
    grid.o 
alphadisk.o: globals.o cgs.o balance.o ss73solution.o rk4settings.o            \
    fileunits.o grid.o odeintegr.o thrfun.o 
dv-alpha.o: settings.o globals.o alphadisk.o fileunits.o results.o 
settings.o: globals.o fileunits.o 
rk4settings.o: globals.o 
relaxation.o: globals.o cgs.o thrfun.o 
modelmag.o: settings.o globals.o cgs.o balance.o ss73solution.o rk4settings.o  \
    odeintegr.o thrfun.o 
globals.o: cgs.o 
