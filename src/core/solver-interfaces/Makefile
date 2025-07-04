# $Id$
#

topdir := ../../..

include $(topdir)/Make.common

# we have to compile each module
subdirs := petsc slepc
ifeq ($(enable_pardiso), yes)
 subdirs += pardiso
endif


include $(topdir)/Make.rules

