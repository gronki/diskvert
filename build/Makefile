#------------------------------------------------------------------------------#
#----------------------------- INSTALLATION PATHS -----------------------------#
#------------------------------------------------------------------------------#

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
libdir = $(exec_prefix)/lib
fmoddir = $(libdir)/finclude/$(fmodsubdir)

#------------------------------------------------------------------------------#
#---------------------------- CFLAGS AND COMPILERS ----------------------------#
#------------------------------------------------------------------------------#

FC := gfortran
FFLAGS := -O2

#------------------------------------------------------------------------------#

FCVER := $(shell $(FC) --version \
	| sed '/^\\s*$$/d' | head -1 | awk '{ print $$1 }')
fmodsubdir := $(notdir $(shell which $(FC)))
LINALG := -llapack -lblas

#------------------------------------------------------------------------------#

ifeq ($(FCVER), GNU)

FFLAGS += -Wall -Wextra -Wno-unused-dummy-argument -pedantic
FFLAGS += -Warray-temporaries -Wrealloc-lhs-all
override FFLAGS += -fimplicit-none -ffree-line-length-none

hasopenblas := $(shell ldconfig -p | grep libopenblas.so > /dev/null \
	&& echo yes || echo no)

ifeq ($(hasopenblas), yes)
$(info OpenBLAS detected! Using instead of BLAS and LAPACK.)
LINALG := -lopenblas
endif

#------------------------------------------------------------------------------#

else ifeq ($(FCVER), ifort)

fmodsubdir := intel
override FFLAGS += -std15 -implicitnone -warn
LINALG := -mkl=sequential

endif

#------------------------------------------------------------------------------#

confortdir := ../libconfort/build
INCLUDE := -I $(confortdir)
LDFLAGS := -L.
LDLIBS += $(LINALG)

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
	install -d $(DESTDIR)$(fmoddir)/diskvert
	install -m 644 ss73solution.mod globals.mod relaxation.mod modelmag.mod \
		alphadisk.mod heatbalance.mod rk4settings.mod slf_*.mod \
		$(DESTDIR)$(fmoddir)/diskvert
	# pkg-config file
	install -d $(DESTDIR)$(libdir)/pkgconfig
	@echo "Name: diskvert" | tee diskvert.pc
	@echo "Description: Set of tools for solving the vertical structure" \
		"of black hole accretion disks." | tee -a diskvert.pc
	@echo "Version: 180706" | tee -a diskvert.pc
	@echo "Libs: -ldiskvert" | tee -a diskvert.pc
	@echo "Cflags: -I$(fmoddir)/diskvert" | tee -a diskvert.pc
	install -m 644 diskvert.pc $(DESTDIR)$(libdir)/pkgconfig
	$(RM) diskvert.pc

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
	$(MAKE) -C $(confortdir) FC=$(FC) libconfort.a

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
