# $Id$
#

topdir := ../
include $(topdir)/Make.common

subdirs := BULK_3D



LDFLAGS += -Wl,-rpath,/usr/pack/gsl-1.6-ma/lib -L/usr/pack/gsl-1.6-ma/lib
LDFLAGS += -lgsl -lgslcblas


.PHONY: all clean

all: target := all
clean: target := clean

all: $(subdirs)

run: all
	for i in $(subdirs);	\
	do			\
	  $(MAKE) -C $$i run	\
	done

clean: $(subdirs)



include $(topdir)/Make.rules

