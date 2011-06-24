# $Id$

topdir := ../../../../..


# we have to find the necessary include files of the module
MODULE_INCLUDES = -I..

MODULE_LDFLAGS = -L$(module_dir)/$(MODULE_NAME)/driftdiffusion$(shared_libext)

include $(topdir)/Make.module

$(foreach lib, $(srcfiles), $(eval $(call MODULE, $(lib))))


include $(topdir)/Make.rules
