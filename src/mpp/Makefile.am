lib_LIBRARIES = libmpp.a

libmpp_a_SOURCES = affinity.c nsclock.c threadloc.c
libmpp_a_SOURCES += mpp.F90 mpp_data.F90 mpp_io.F90 mpp_memutils.F90 mpp_pset.F90 mpp_domains.F90 \
	mpp_efp.F90 mpp_parameter.F90 mpp_utilities.F90

EXTRA_DIST = fms include

local_INCLUDE = -I . -I $(top_srcdir)/src/mpp/include/ -I $(top_srcdir)/src/mpp/fms/
libmpp_a_INCLUDE  = $(local_INCLUDE)
libmpp_a_FCFLAGS  = $(local_INCLUDE)
libmpp_a_FCFLAGS += -DOVERLOAD_R4 -Duse_libMPI

clean-local:
	-rm -f *.mod

fortran-local: all
nuopc-local: all

libmpp_a-affinity.$(OBJEXT) :
libmpp_a-mpp_parameter.$(OBJEXT) : fms/fms_platform.h fms/file_version.h
libmpp_a-mpp_data.$(OBJEXT)      : fms/fms_platform.h fms/file_version.h include/mpp_data_sma.inc include/mpp_data_mpi.inc include/mpp_data_nocomm.inc
libmpp_a-mpp_data.$(OBJEXT)    : libmpp_a-mpp_parameter.$(OBJEXT)
libmpp_a-mpp_domains.$(OBJEXT) : fms/fms_platform.h fms/file_version.h include/mpp_define_nest_domains.inc include/mpp_domains_util.inc \
	include/mpp_domains_comm.inc include/mpp_domains_define.inc include/mpp_domains_misc.inc include/mpp_domains_reduce.inc
libmpp_a-mpp_domains.$(OBJEXT) : libmpp_a-mpp_parameter.$(OBJEXT) libmpp_a-mpp_data.$(OBJEXT) libmpp_a-mpp.$(OBJEXT) libmpp_a-mpp_memutils.$(OBJEXT) \
	libmpp_a-mpp_pset.$(OBJEXT) libmpp_a-mpp_efp.$(OBJEXT)
libmpp_a-mpp_efp.$(OBJEXT)     : fms/fms_platform.h 
libmpp_a-mpp_efp.$(OBJEXT)     : libmpp_a-mpp.$(OBJEXT)
libmpp_a-mpp.$(OBJEXT)         : fms/fms_platform.h fms/file_version.h include/system_clock.h include/mpp_util.inc include/mpp_comm.inc
libmpp_a-mpp.$(OBJEXT)         : libmpp_a-mpp_parameter.$(OBJEXT) libmpp_a-mpp_data.$(OBJEXT) fms/fms_platform.h fms/file_version.h
libmpp_a-mpp_io.$(OBJEXT)      : fms/fms_platform.h fms/file_version.h include/mpp_io_util.inc include/mpp_io_misc.inc include/mpp_io_connect.inc \
	include/mpp_io_read.inc include/mpp_io_write.inc
libmpp_a-mpp_io.$(OBJEXT)      : libmpp_a-mpp_parameter.$(OBJEXT) libmpp_a-mpp.$(OBJEXT) libmpp_a-mpp_domains.$(OBJEXT) 
libmpp_a-mpp_memutils.$(OBJEXT): libmpp_a-mpp.$(OBJEXT)
libmpp_a-mpp_pset.$(OBJEXT)    : fms/fms_platform.h
libmpp_a-mpp_pset.$(OBJEXT)    : libmpp_a-mpp.$(OBJEXT)
libmpp_a-mpp_utilities.$(OBJEXT) : fms/fms_platform.h
libmpp_a-mpp_utilities.$(OBJEXT) : libmpp_a-mpp.$(OBJEXT)
