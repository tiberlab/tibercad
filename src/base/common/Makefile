# $Id$

topdir := ../../..

include $(topdir)/Make.common

all: svnrevision


.PHONY: svnrevision
svnrevision:
	@if test ! -e $@.h; then echo "#define SVNREVISION $(SVNREVISION)" > $@.h; fi; 
	@if test `awk '{print $$3}' $@.h` -ne $(SVNREVISION); then echo "#define SVNREVISION $(SVNREVISION)" > $@.h; fi



include $(topdir)/Make.rules
