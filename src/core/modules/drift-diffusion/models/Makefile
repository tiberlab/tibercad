# $Id$

topdir := ../../../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
ifeq ($(build_modules), yes)
 objects	:= $(patsubst %.C, %.so, $(srcfiles))
else
 objects	:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
endif

CXXFLAGS += $(MODULE_CXXFLAGS)


.PHONY: all clean

all: $(objects)


%.so : %.$(obj-suffix) .%.d
	@soname=lib`$(CPP) $(CPPFLAGS) -dM $(@:.so=.C) \
		| awk '/TIBER_MODULE_NAME/ { module = $$3; } \
		  	END { if (module != "") printf("%s_", module) }'`; \
				soname=$$soname`grep TIBER_MODULE $(@:.so=.C) | \
					sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`.so ; \
		echo "Linking module $$soname ..."; \
		$(CXX) -shared  $(CPPFLAGS) $(CXXFLAGS) \
		-Wl,-soname,$$soname -o $@ $<; \
		cp $@ $(topdir)/lib/tibermodels/$$soname



clean:
	@rm -f *.so




include $(topdir)/Make.rules
