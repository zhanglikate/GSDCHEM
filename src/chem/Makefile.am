SUBDIRS = shr plume vash dust aero dep seas opt gocart .

lib_LIBRARIES = libgsdchem.a

include $(ESMFMKFILE)

libgsdchem_a_LIBADD  = shr/*.$(OBJEXT) plume/*.$(OBJEXT) vash/*.$(OBJEXT) dust/*.$(OBJEXT) dep/*.$(OBJEXT) aero/*.$(OBJEXT) seas/*.$(OBJEXT) opt/*.$(OBJEXT) gocart/*.$(OBJEXT)

#local_INCLUDE      = -I . -I $(top_builddir)/src/mpp -I $(top_builddir)/src/util/ -I shr/ -I gocart -I plume -I vash -I dust -I dep -I aero -I seas -I opt
local_INCLUDE      = -I . -I shr/ -I gocart -I plume -I vash -I dust -I dep -I aero -I seas -I opt
libgsdchem_a_INCLUDE  = $(local_INCLUDE)
libgsdchem_a_FCFLAGS  = $(local_INCLUDE)
libgsdchem_a_FCFLAGS += $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS)

libgsdchem_a_SOURCES = chemCap.F90 chem_methods.F90

clean-local:
	-rm -f *.mod

libgsdchem_a-chemCap.$(OBJEXT): libgsdchem_a-chem_methods.$(OBJEXT)
libgsdchem_a-chem_methods.$(OBJEXT): $(builddir)/shr/libshr.a
