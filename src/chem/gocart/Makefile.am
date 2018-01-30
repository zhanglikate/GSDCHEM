noinst_LIBRARIES = libgocart.a

libgocart_a_SOURCES = gocart_mod.F90 
libgocart_a_SOURCES += gocart_settling_mod.F90
libgocart_a_SOURCES += gocart_aerosols_mod.F90
libgocart_a_SOURCES += gocart_dmsemis_mod.F90
libgocart_a_SOURCES += gocart_chem_mod.F90

local_INCLUDE =
local_INCLUDE += -I .
#local_INCLUDE += -I $(top_builddir)/src/util
local_INCLUDE += -I ../shr/ -I ../plume/ -I ../vash/ -I ../dust/ -I ../dep/ -I ../seas/ -I ../opt/ -I ../aero/

libgocart_a_INCLUDE = $(local_INCLUDE)
libgocart_a_FCFLAGS = $(local_INCLUDE)

#libgocart_a_DEPENDENCIES = $(top_builddir)/src/util/libtime.a
libgocart_a_DEPENDENCIES = ../shr/libshr.a ../plume/libplume.a ../vash/libvash.a ../dust/libdust.a ../dep/libdep.a ../seas/libseas.a ../opt/libopt.a ../aero/libaero.a

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all

libgocart_a-gocart_mod.$(OBJEXT) : libgocart_a-gocart_settling_mod.$(OBJEXT) libgocart_a-gocart_aerosols_mod.$(OBJEXT) \
	libgocart_a-gocart_dmsemis_mod.$(OBJEXT) libgocart_a-gocart_chem_mod.$(OBJEXT) 
