noinst_LIBRARIES = libdust.a

libdust_a_SOURCES = dust_afwa_mod.F90 dust_data_mod.F90 dust_gocart_mod.F90 dust_mod.F90

#local_INCLUDE = -I ../shr/
#libdust_a_INCLUDE = $(local_INCLUDE)
libdust_a_FCFLAGS = -I ../shr/

libdust_a_DEPENDENCIES = ../shr/libshr.a

clean-local:
	-rm -f *.mod

libdust_a-dust_afwa_mod.$(OBJEXT): libdust_a-dust_data_mod.$(OBJEXT)
libdust_a-dust_gocart_mod.$(OBJEXT): libdust_a-dust_data_mod.$(OBJEXT)
libdust_a-dust_mod.$(OBJEXT): libdust_a-dust_data_mod.$(OBJEXT) libdust_a-dust_gocart_mod.$(OBJEXT) libdust_a-dust_afwa_mod.$(OBJEXT)
