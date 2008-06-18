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



.PHONY: all clean

all: $(objects)

test:
	libtool --mode=compile g++ $(CPPFLAGS) $(CXXFLAGS) -c ConstantMobility.C
	libtool --mode=link g++ -module $(LDFLAGS) -o libdd_mob_constant.la ConstantMobility.lo -rpath $(topdir)/lib/tibermodels -dlopen libunstrained.la


%.so : %.$(obj-suffix)
	@soname=lib`$(CPP) $(CPPFLAGS) $(CXXFLAGS) -dM $(@:.so=.C) \
		| awk '/TIBER_MODULE_NAME/ { module = $$3; } \
		  	END { if (module != "") printf("%s_", module) }'`; \
				soname=$$soname`grep TIBER_MODULE $(@:.so=.C) | \
					sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`.so ; \
		echo "Linking module $$soname ..."; \
		$(CXX) -shared $(CXXFLAGS) -Wl,-soname,$$soname -o $@ $<; \
		cp $@ $(topdir)/lib/tibermodels/$$soname



clean:
	@rm -f *.so *~ *.d .*.d .*.d.*



ifneq ($(MAKECMDGOALS),clean)
-include $(patsubst %,.%,$(srcfiles:.C=.d))
endif

include $(topdir)/Make.rules
