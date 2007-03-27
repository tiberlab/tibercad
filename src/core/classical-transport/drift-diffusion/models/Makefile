# $Id$

topdir := ../../../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))



.PHONY: all clean

all: $(objects)

test:
	@$(CPP) $(CPPFLAGS) -DTIBER_MODULE -dM ConstantMobility.C \
		| awk '/TIBER_MODULENAME/ { module = $$3; } \
		  /TIBER_MODELNAME/ { model = $$3; } \
		  END { if (module == "") printf("%s", model); \
		        else printf("%s_%s", module, model); }'


clean:
	@rm -f *.o *~ *.d .*.d .*.d.*



ifneq ($(MAKECMDGOALS),clean)
-include $(patsubst %,.%,$(srcfiles:.C=.d))
endif

include $(topdir)/Make.rules
