noinst_LIBRARIES = libshr.a

libshr_a_SOURCES = chem_comm_mod.F90 chem_config_mod.F90 chem_const_mod.F90 chem_domain_mod.F90 chem_io_mod.F90 chem_shr_mod.F90 chem_species_mod.F90 chem_state_mod.F90 chem_types_mod.F90 chem_vars_mod.F90 chem_io_mod.F90 units.F90

local_INCLUDE  =
local_INCLUDE += -I .
local_INCLUDE += -I $(top_builddir)/src/mpp
local_INCLUDE += -I $(top_builddir)/src/util

libshr_a_INCLUDE = $(local_INCLUDE)
libshr_a_FCFLAGS = $(local_INCLUDE)

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all

libshr_a-chem_config_mod.$(OBJEXT) : libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_const_mod.$(OBJEXT)  : libshr_a-chem_types_mod.$(OBJEXT)
libshr_a-chem_io_mod.$(OBJEXT)     : libshr_a-chem_domain_mod.$(OBJEXT) libshr_a-chem_config_mod.$(OBJEXT)
libshr_a-chem_species_mod.$(OBJEXT): libshr_a-chem_state_mod.$(OBJEXT) libshr_a-chem_config_mod.$(OBJEXT)
libshr_a-chem_shr_mod.$(OBJEXT)    : libshr_a-chem_const_mod.$(OBJEXT) libshr_a-chem_state_mod.$(OBJEXT) \
	libshr_a-chem_vars_mod.$(OBJEXT) libshr_a-chem_config_mod.$(OBJEXT) libshr_a-chem_domain_mod.$(OBJEXT) \
	libshr_a-chem_species_mod.$(OBJEXT) libshr_a-chem_io_mod.$(OBJEXT) 
