SUBDIRS = src
dist_doc_DATA = README
EXTRA_DIST = test

install-data-hook: $(datarootdir)/gsdchem.mk

$(datarootdir)/gsdchem.mk:
	@echo "# ESMF self-describing build dependency makefile fragment" > $@
	@echo >> $@
	@echo "ESMF_DEP_FRONT     = CHM"                    >> $@
	@echo "ESMF_DEP_INCPATH   = $(abs_top_builddir)/src/chem/" >> $@
	@echo "ESMF_DEP_CMPL_OBJS = "                       >> $@
	@echo "ESMF_DEP_LINK_OBJS = $(libdir)/libgsdchem.a" >> $@
	@echo "ESMF_DEP_SHRD_PATH = "                   >> $@
	@echo "ESMF_DEP_SHRD_LIBS = mpp"                >> $@
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" $@
	@echo
