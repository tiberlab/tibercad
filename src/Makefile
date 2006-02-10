# $Id$
#

topdir := ../
include $(topdir)/Make.common

subdirs := utils/Tensor base/io base/common core/drift-diffusion core/macrostrain



.PHONY: all clean

all: target := all
clean: target := clean

all: $(subdirs)

clean: $(subdirs)



include $(topdir)/Make.rules

