# $Id$
#

topdir := ../..

include $(topdir)/Make.common

# we have to compile each module
subdirs := Tensor SLEPC 
subdirs += license

ifeq ($(enable_hetero), yes)
 subdirs += HETERO
endif


include $(topdir)/Make.rules

