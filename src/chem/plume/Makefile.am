noinst_LIBRARIES = libplume.a

libplume_a_SOURCES = plume_rise_mod.F90 plume_zero_mod.F90 plume_scalar_mod.F90

#local_INCLUDE      = -I ../shr/ -I $(top_builddir)/src/mpp
local_INCLUDE      = -I ../shr/
libplume_a_INCLUDE = $(local_INCLUDE)
libplume_a_FCFLAGS = $(local_INCLUDE)

#libplume_a_DEPENDENCIES = ../shr/libshr.a $(top_builddir)/src/mpp/libmpp.a
libplume_a_DEPENDENCIES = ../shr/libshr.a

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all

libplume_a-plume_rise_mod.$(OBJEXT)   : libplume_a-plume_zero_mod.$(OBJEXT) libplume_a-plume_scalar_mod.$(OBJEXT)
libplume_a-plume_scalar_mod.$(OBJEXT) : libplume_a-plume_zero_mod.$(OBJEXT)
