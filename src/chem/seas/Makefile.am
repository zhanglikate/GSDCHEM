noinst_LIBRARIES = libseas.a

libseas_a_SOURCES = seas_data_mod.F90 seas_mod.F90

local_INCLUDE     = -I ../shr/ 
libseas_a_INCLUDE = $(local_INCLUDE)
libseas_a_FCFLAGS = $(local_INCLUDE)

libseas_a_DEPENDENCIES = ../shr/libshr.a

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all

libseas_a-seas_mod.$(OBJEXT)   : libseas_a-seas_data_mod.$(OBJEXT)
