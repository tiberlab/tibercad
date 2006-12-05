# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))



.PHONY: all clean


all: $(objects)


clean:
	@rm -f *.o *~ *.d .*.d .*.d.*



ifneq ($(MAKECMDGOALS),clean)
-include $(patsubst %,.%,$(srcfiles:.C=.d))
endif

include $(topdir)/Make.rules
