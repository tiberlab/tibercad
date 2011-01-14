# $Id$

topdir := ../../..

include $(topdir)/Make.common

all: svnrevision.h


.PHONY: svnrevision.h
svnrevision.h:
	@if test ! -e $@; then echo "#define SVNREVISION $(SVNREVISION)" > $@; fi; 
	@if test `awk '{print $$3}' $@` -ne $(SVNREVISION); then echo "#define SVNREVISION $(SVNREVISION)" > $@; fi



include $(topdir)/Make.rules
