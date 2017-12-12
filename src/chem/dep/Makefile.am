noinst_LIBRARIES = libdep.a

libdep_a_SOURCES = dep_mod.F90 dep_wet_ls_mod.F90 dep_dry_gocart_mod.F90 dep_dry_mod.F90 dep_simple_mod.F90 \
	dep_vertmx_mod.F90

local_INCLUDE    = -I . -I ../shr/ -I ../aero/
libdep_a_INCLUDE = $(local_INCLUDE)
libdep_a_FCFLAGS = $(local_INCLUDE)

libdep_a_DEPENDENCIES = ../shr/libshr.a ../aero/libaero.a 

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all

libdep_a-dep_mod.$(OBJEXT)     : libdep_a-dep_wet_ls_mod.$(OBJEXT) libdep_a-dep_dry_mod.$(OBJEXT)
libdep_a-dep_dry_mod.$(OBJEXT) : libdep_a-dep_dry_gocart_mod.$(OBJEXT) libdep_a-dep_simple_mod.$(OBJEXT) libdep_a-dep_vertmx_mod.$(OBJEXT)
