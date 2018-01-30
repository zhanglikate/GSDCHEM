noinst_LIBRARIES = libaero.a

libaero_a_SOURCES = aero_soa_vbs_mod.F90 aero_soa_vbs_data_mod.F90

local_INCLUDE     = -I . -I ../shr/ -I ../dust/
libaero_a_INCLUDE = $(local_INCLUDE)
libaero_a_FCFLAGS = $(local_INCLUDE)

libaero_a_DEPENDENCIES = ../shr/libshr.a ../dust/libdust.a

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all

libaero_a-aero_soa_vbs_mod.$(OBJEXT)   : libaero_a-aero_soa_vbs_data_mod.$(OBJEXT)
