#------------------------------------------------------------------------------#
#----------------------------- INSTALLATION PATHS -----------------------------#
#------------------------------------------------------------------------------#

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
libdir = $(exec_prefix)/lib

ifeq ($(fmoddir),)
fmoddir_ = $(libdir)/diskvert/modules
else
fmoddir_ = $(fmoddir)/diskvert
endif

#------------------------------------------------------------------------------#
#---------------------------- CFLAGS AND COMPILERS ----------------------------#
#------------------------------------------------------------------------------#

FC := gfortran
FFLAGS = -O2 -funsafe-math-optimizations

override FFLAGS += -Wall -Wno-unused-dummy-argument \
		-Warray-temporaries -Wrealloc-lhs-all \
		-pedantic -fimplicit-none -ffree-line-length-none

#------------------------------------------------------------------------------#

confortdir := ../libconfort/build
INCLUDE := -I $(confortdir)
LDFLAGS := -L.
LDLIBS += -lopenblas

#------------------------------------------------------------------------------#
#-------------------------- SOURCES AND OBJECT FILES --------------------------#
#------------------------------------------------------------------------------#

# search for sources in these directories
VPATH := ../src:../src/util:../src/prog:../src/math
programs = dv-alpha dv-alpha-rx dv-mag dv-mag-rx dv-rad1

objects_math := cgs.o findzer.o integrderiv.o interpol.o odeintegr.o thrfun.o
objects_util := fileunits.o results.o settings.o
objects_disk := globals.o balance.o grid.o ss73solution.o
objects_lib = $(objects_math) $(objects_disk) alphadisk.o modelmag.o rk4settings.o relaxation.o
objects_comm = $(objects_math) $(objects_disk) $(objects_util)

all: $(programs) libdiskvert.so

#------------------------------------------------------------------------------#
#------------------------------ INSTALLATION ----------------------------------#
#------------------------------------------------------------------------------#

install: all
	install -d $(DESTDIR)$(bindir)
	install $(programs) $(DESTDIR)$(bindir)
	ln -sf $(bindir)/dv-mag-rx $(DESTDIR)$(bindir)/diskvert
	# library file
	install -d $(DESTDIR)$(libdir)
	install libdiskvert.so $(DESTDIR)$(libdir)
	# Fortran modules
	install -d $(DESTDIR)$(fmoddir_)
	install -m 644 ss73solution.mod globals.mod relaxation.mod modelmag.mod \
		alphadisk.mod heatbalance.mod rk4settings.mod slf_*.mod \
		$(DESTDIR)$(fmoddir_)
	# pkg-config file
	install -d $(DESTDIR)$(libdir)/pkgconfig
	@echo "Name: diskvert" | tee diskvert.pc
	@echo "Description: Set of tools for solving the vertical structure" \
		"of black hole accretion disks." | tee -a diskvert.pc
	@echo "Version: 181031" | tee -a diskvert.pc
	@echo "Libs: -ldiskvert" | tee -a diskvert.pc
	@echo "Cflags: -I$(fmoddir_)" | tee -a diskvert.pc
	install -m 644 diskvert.pc $(DESTDIR)$(libdir)/pkgconfig
	$(RM) diskvert.pc

install-home: override prefix = $(HOME)/.local
install-home: install

#------------------------------------------------------------------------------#
#-------------------------------- BUILD RULES ---------------------------------#
#------------------------------------------------------------------------------#

%.o: %.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c $< -o $@

include deps.inc

relaxation.o: mrxcoeff.fi mrxdims.fi mrxhash.fi mrxptrs.fi
settings.o: $(confortdir)/libconfort.a

#------------------------------------------------------------------------------#

$(programs): LDFLAGS += -fwhole-program

dv-alpha-rx: $(objects_comm) relaxation.o dv-alpha-rx.f90 $(confortdir)/libconfort.a
	$(FC) $(INCLUDE) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

dv-mag-rx: $(objects_comm) relaxation.o dv-mag-rx.f90 $(confortdir)/libconfort.a
	$(FC) $(INCLUDE) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

dv-alpha: $(objects_comm) alphadisk.o rk4settings.o dv-alpha.f90 $(confortdir)/libconfort.a
	$(FC) $(INCLUDE) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

dv-mag: $(objects_comm) modelmag.o rk4settings.o dv-mag.f90 $(confortdir)/libconfort.a
	$(FC) $(INCLUDE) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

dv-rad1: $(objects_comm) alphadisk.o rk4settings.o dv-rad1.f90 $(confortdir)/libconfort.a
	$(FC) $(INCLUDE) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

$(objects_lib): override FFLAGS += -fpic

libdiskvert.so: $(objects_lib)
	$(FC) $(FFLAGS) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

$(confortdir)/libconfort.a:
	$(MAKE) -C $(confortdir) FC=$(FC) CFLAGS='-O2' FFLAGS='-O2' libconfort.a

#------------------------------------------------------------------------------#
#---------------------------------- CLEANUP -----------------------------------#
#------------------------------------------------------------------------------#

clean:
	$(RM) *.mod *.smod *.a *.o $(programs) *.so
	$(MAKE) -C $(confortdir) clean

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

.PHONY: all install clean

#------------------------------------------------------------------------------#
