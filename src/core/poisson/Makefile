# $Id: Makefile 163 2007-01-12 12:37:28Z maufder $

topdir := ../../..

include $(topdir)/Make.common


#ifeq ($(build_modules), yes)
# CXXFLAGS += $(MODULE_CXXFLAGS)
# target := shared
#endif



shared: $(objects)
	@soname=lib`grep TIBER_MODULE Poisson.C | \
		sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`.so ; \
		echo "Linking module $$soname ..."; \
		$(CXX) -shared $(CXXFLAGS) -Wl,-soname,$$soname \
			-o $(libdir)/tibermodels/$$soname $(objects)




include $(topdir)/Make.rules
