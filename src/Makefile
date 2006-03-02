# $Id$
#

topdir := ../
include $(topdir)/Make.common

# we have to compile each module
subdirs := $(tibermodules)


.PHONY: all clean

all: target := all
clean: target := clean

all: $(subdirs)

clean: $(subdirs)



include $(topdir)/Make.rules

