noinst_LIBRARIES = libgocart.a

libgocart_a_SOURCES  = gocart_mod.F90
libgocart_a_SOURCES += gocart_aerosols_mod.F90
libgocart_a_SOURCES += gocart_chem_mod.F90
libgocart_a_SOURCES += gocart_dmsemis_mod.F90
libgocart_a_SOURCES += gocart_model_mod.F90
libgocart_a_SOURCES += gocart_prep_mod.F90
libgocart_a_SOURCES += gocart_settling_mod.F90

libgocart_a_FCFLAGS = -I . -I ../shr/ -I ../plume/ -I ../vash/ -I ../dust/ -I ../dep/ -I ../seas/ -I ../opt/

libgocart_a_DEPENDENCIES = ../shr/libshr.a ../plume/libplume.a ../vash/libvash.a ../dep/libdep.a ../dust/libdust.a ../seas/libseas.a ../opt/libopt.a

clean-local:
	-rm -f *.mod

libgocart_a-gocart_mod.$(OBJEXT) : libgocart_a-gocart_aerosols_mod.$(OBJEXT) libgocart_a-gocart_dmsemis_mod.$(OBJEXT) \
	libgocart_a-gocart_chem_mod.$(OBJEXT) libgocart_a-gocart_prep_mod.$(OBJEXT) libgocart_a-gocart_settling_mod.$(OBJEXT)
libgocart_a-gocart_model_mod.$(OBJEXT) : libgocart_a-gocart_mod.$(OBJEXT)
