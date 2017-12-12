noinst_LIBRARIES = libvash.a

libvash_a_SOURCES = vash_settling_mod.F90

local_INCLUDE     = -I ../shr/ 
libvash_a_INCLUDE = $(local_INCLUDE)
libvash_a_FCFLAGS = $(local_INCLUDE)

libvash_a_DEPENDENCIES = ../shr/libshr.a 

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all
