noinst_LIBRARIES = libseas.a

libseas_a_SOURCES = seas_data_mod.F90 seas_mod.F90 seas_ngac_mod.F90

libseas_a_FCFLAGS = -I ../shr/

libseas_a_DEPENDENCIES = ../shr/libshr.a

clean-local:
	-rm -f *.mod

libseas_a-seas_mod.$(OBJEXT) : libseas_a-seas_data_mod.$(OBJEXT) libseas_a-seas_ngac_mod.$(OBJEXT)
