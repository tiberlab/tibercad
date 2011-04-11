# $Id$
#

topdir := ../..

include $(topdir)/Make.common

# we have to compile each module
subdirs := SLEPC Tensor license
ifeq ($(enable_dftb), yes)
 subdirs += DFTB
endif
ifeq ($(enable_uptight), yes)
 subdirs += UPTIGHT
endif
ifeq ($(enable_hetero), yes)
 subdirs += HETERO
endif


include $(topdir)/Make.rules

