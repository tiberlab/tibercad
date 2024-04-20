# $Id$
#

topdir := ../..

include $(topdir)/Make.common

# we have to compile each module
subdirs := Tensor SLEPC 
ifeq ($(enable_license_check), yes)
subdirs += license
endif
ifeq ($(enable_hetero), yes)
 subdirs += HETERO
endif


include $(topdir)/Make.rules

