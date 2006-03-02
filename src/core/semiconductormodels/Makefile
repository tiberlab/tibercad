# $Id$

topdir := ../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))



.PHONY: all clean


all: $(objects)

#-include $(srcfiles:.C=.d)


clean:
	@rm -f *.o *~ *.d



