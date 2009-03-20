# $Id$

topdir := ../../../..

# subdirectories
subdirs		:= 


include $(topdir)/Make.common

ifeq ($(build_modules), yes)
 target := shared
endif



shared: $(objects)
	@soname=lib`grep TIBER_MODULE DSSC.C | \
		sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`$(shared_libext) ; \
		echo "Linking module $$soname ..."; \
		$(CXX) -shared $(CXXFLAGS) -Wl,-soname,$$soname \
			-o $(libdir)/tibermodels/$$soname $(objects)


	

include $(topdir)/Make.rules
