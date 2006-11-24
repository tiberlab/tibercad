# $Id$
#

topdir := ../
include $(topdir)/Make.common

subdirs := BULK MESFET



.PHONY: all clean

# sets the target for the subdirs
all: target := all
clean: target := clean
distclean: target := distclean

all: $(subdirs)

run: all
	for i in $(subdirs);	\
	do			\
	  $(MAKE) -C $$i run	\
	done

clean: $(subdirs)

distclean: $(subdirs)



include $(topdir)/Make.rules

