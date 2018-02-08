noinst_LIBRARIES = libplume.a

libplume_a_SOURCES = plume_rise_mod.F90 plume_zero_mod.F90 plume_scalar_mod.F90

libplume_a_FCFLAGS = -I ../shr/

libplume_a_DEPENDENCIES = ../shr/libshr.a

clean-local:
	-rm -f *.mod

libplume_a-plume_rise_mod.$(OBJEXT)   : libplume_a-plume_zero_mod.$(OBJEXT) libplume_a-plume_scalar_mod.$(OBJEXT)
libplume_a-plume_scalar_mod.$(OBJEXT) : libplume_a-plume_zero_mod.$(OBJEXT)
