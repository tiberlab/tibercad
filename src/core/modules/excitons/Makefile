# $Id$

topdir := ../../../..

include $(topdir)/Make.common


#
# source files
mainsrcfiles	:= ExcitonTransport.C ExcitonProperties.C
modelsrcfiles	:= ExcitonModel.C

#
# object files
mainobjects	:= $(patsubst %.C, %$(obj-suffix), $(mainsrcfiles))
modelobjects	:= $(patsubst %.C, %$(obj-suffix), $(modelsrcfiles))

ifeq ($(build_modules), yes)
 target := shared $(patsubst %.C, %$(shared_libext), $(modelsrcfiles))
else
 target := $(mainobjects) $(modelobjects)
endif



shared: $(mainobjects)
	@soname=lib`grep TIBER_MODULE ExcitonTransport.C | \
		sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`$(shared_libext) ; \
		echo "Linking module $$soname ..."; \
		$(CXX) -shared $(CXXFLAGS) -Wl,-soname,$$soname \
			-o $(libdir)/tibermodels/$$soname $(mainobjects)


%$(shared_libext) : %.C
	@soname=lib`$(CPP) $(CPPFLAGS) -dM $(@:$(shared_libext)=.C) \
		| awk '/TIBER_MODULE_NAME/ { module = $$3; } \
		  	END { if (module != "") printf("%s_", module) }'`; \
				soname=$$soname`grep TIBER_MODULE $(@:$(shared_libext)=.C) | \
					sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`$(shared_libext) ; \
		echo "Linking module $$soname ..."; \
		$(CXX) -shared $(CPPFLAGS) $(CXXFLAGS) -Wl,-soname,$$soname -o $@ $<; \
		cp $@ $(libdir)/tibermodels/$$soname



clean:
	@rm -f *$(shared_libext)




include $(topdir)/Make.rules
