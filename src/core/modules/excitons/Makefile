# $Id$

topdir := ../../../..

include $(topdir)/Make.common


#
# source files
mainsrcfiles	:= ExcitonTransport.C ExcitonProperties.C
modelsrcfiles	:= ExcitonModel.C

#
# object files
mainobjects	:= $(patsubst %.C, %.$(obj-suffix), $(mainsrcfiles))
modelobjects	:= $(patsubst %.C, %.$(obj-suffix), $(modelsrcfiles))

ifeq ($(build_modules), yes)
 target := shared $(patsubst %.C, %.so, $(modelsrcfiles))
else
 target := $(mainobjects) $(modelobjects)
endif



.PHONY: all
all: $(target)


shared: $(mainobjects)
	@soname=lib`grep TIBER_MODULE ExcitonTransport.C | \
		sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`.so ; \
		echo "Linking module $$soname ..."; \
		$(CXX) -shared $(CXXFLAGS) -Wl,-soname,$$soname \
			-o $(topdir)/lib/tibermodels/$$soname $(mainobjects)


%.so : %.C
	@soname=lib`$(CPP) $(CPPFLAGS) -dM $(@:.so=.C) \
		| awk '/TIBER_MODULE_NAME/ { module = $$3; } \
		  	END { if (module != "") printf("%s_", module) }'`; \
				soname=$$soname`grep TIBER_MODULE $(@:.so=.C) | \
					sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`.so ; \
		echo "Linking module $$soname ..."; \
		$(CXX) -shared $(CPPFLAGS) $(CXXFLAGS) -Wl,-soname,$$soname -o $@ $<; \
		cp $@ $(topdir)/lib/tibermodels/$$soname



clean:
	@rm -f *.so




include $(topdir)/Make.rules
