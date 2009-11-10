# $Id$

topdir := ../../..

include $(topdir)/Make.common



TiberCad$(obj-suffix): svnrevision.h

.PHONY: svnrevision.h
svnrevision.h:
	@echo "#define SVNREVISION $(SVNREVISION)" > svnrevision.h



include $(topdir)/Make.rules
