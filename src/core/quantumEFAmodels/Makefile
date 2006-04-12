# $Id: Makefile 1 2006-02-10 09:27:03Z maufder $

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



