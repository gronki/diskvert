#----------------------------- INSTALLATION PATHS -----------------------------#

prefix := $(shell if test -w "/usr/local"; then echo "/usr/local"; \
        else echo "$(HOME)/.local"; fi)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
libdir = $(exec_prefix)/lib

ifeq ($(fmoddir),)
  fmoddir_ = $(libdir)/diskvert/modules
else
  fmoddir_ = $(fmoddir)/diskvert
endif

#---------------------------- CFLAGS AND COMPILERS ----------------------------#

FC := gfortran
ARCHFLAGS ?= -march=native
FFLAGS := -O3 -funsafe-math-optimizations $(ARCHFLAGS)
FFLAGS_EXTRA := -fimplicit-none -ffree-line-length-none -std=f2018 \
  -pedantic -Wall -Wno-unused-dummy-argument -Warray-temporaries

#------------------------------------------------------------------------------#

confortdir := ../libconfort/build
INCLUDE += -I $(confortdir)
LDLIBS += -lopenblas

#-------------------------- SOURCES AND OBJECT FILES --------------------------#

# search for sources in these directories
VPATH := ../src:../src/util:../src/prog:../src/math
programs = dv-alpha dv-mag diskvert dv-rad1 dv-ss73

objects_math := cgs.o findzer.o integrderiv.o interpol.o odeintegr.o thrfun.o
objects_util := fileunits.o results.o settings.o
objects_disk := globals.o balance.o grid.o ss73solution.o
objects_lib = $(objects_math) $(objects_disk) alphadisk.o modelmag.o rk4settings.o relaxation.o
objects_comm = $(objects_math) $(objects_disk) $(objects_util)

all: $(programs) libdiskvert.so

#------------------------------ INSTALLATION ----------------------------------#

install: all
	install -d $(DESTDIR)$(bindir)
	install $(programs) $(DESTDIR)$(bindir)
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
	@echo "Version: 200312" | tee -a diskvert.pc
	@echo "Libs: -ldiskvert" | tee -a diskvert.pc
	@echo "Cflags: -I$(fmoddir_)" | tee -a diskvert.pc
	install -m 644 diskvert.pc $(DESTDIR)$(libdir)/pkgconfig
	$(RM) diskvert.pc

#-------------------------------- BUILD RULES ---------------------------------#

%.o: %.f90
	$(FC) $(INCLUDE) $(FFLAGS) -fpic $(FFLAGS_EXTRA) -c $< -o $@

include deps.inc
settings.o $(programs:=.o): $(confortdir)/libconfort.a

#------------------------------------------------------------------------------#

$(programs): LDFLAGS += -fwhole-program
$(programs): LDLIBS += $(confortdir)/libconfort.a

$(programs):
	$(FC) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

libdiskvert.so: $(objects_lib)
	$(FC) $(FFLAGS) $(LDFLAGS) -shared $^ $(LDLIBS) -o $@

#------------------------------------------------------------------------------#

$(confortdir)/libconfort.a:
	$(MAKE) -C $(confortdir) FC="$(FC)" libconfort.a

#---------------------------------- CLEANUP -----------------------------------#

clean:
	$(RM) *.mod *.smod *.a *.o $(programs) *.so
	$(MAKE) -C $(confortdir) clean

#------------------------------------------------------------------------------#

.PHONY: all install clean confort

#------------------------------------------------------------------------------#
