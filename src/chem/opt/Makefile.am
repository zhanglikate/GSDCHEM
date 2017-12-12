noinst_LIBRARIES = libopt.a

libopt_a_SOURCES = opt_averaging_mod.F90 opt_data_mod.F90 opt_gocart_mod.F90 opt_driver_mod.F90 opt_mod.F90
libopt_a_SOURCES += opt_peg_util_mod.F90
libopt_a_SOURCES += opt_aer_mod.F90 opt_aer_out_mod.F90 opt_aer_ra_mod.F90

local_INCLUDE    = -I . -I $(top_builddir)/src/mpp -I ../shr/ -I ../aero/ -I ../dust/ -I ../seas/
libopt_a_INCLUDE = $(local_INCLUDE)
libopt_a_FCFLAGS = $(local_INCLUDE)

libopt_a_DEPENDENCIES = ../shr/libshr.a ../aero/libaero.a ../dust/libdust.a ../seas/libseas.a $(top_builddir)/src/mpp

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all

libopt_a-opt_mod.$(OBJEXT)   : libopt_a-opt_data_mod.$(OBJEXT) libopt_a-opt_averaging_mod.$(OBJEXT) libopt_a-opt_driver_mod.$(OBJEXT) libopt_a-opt_gocart_mod.$(OBJEXT) libopt_a-opt_aer_mod.$(OBJEXT)
libopt_a-opt_driver_mod.$(OBJEXT)   : libopt_a-opt_data_mod.$(OBJEXT) libopt_a-opt_averaging_mod.$(OBJEXT)
libopt_a-opt_averaging_mod.$(OBJEXT) : libopt_a-opt_data_mod.$(OBJEXT) libopt_a-opt_peg_util_mod.$(OBJEXT)
libopt_a-opt_aer_mod.$(OBJEXT) : libopt_a-opt_aer_out_mod.$(OBJEXT) libopt_a-opt_aer_ra_mod.$(OBJEXT)
