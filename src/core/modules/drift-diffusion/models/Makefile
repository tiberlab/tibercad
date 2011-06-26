# $Id$

topdir := ../../../../..


# we have to find the necessary include files of the module
MODULE_INCLUDES = -I..

include $(topdir)/Make.module

MODULE_LIB := $(module_dir)/$(MODULE_NAME)/driftdiffusion$(shared_libext)

$(foreach lib, $(srcfiles), $(eval $(call MODULE, $(lib), $(MODULE_LIB))))


include $(topdir)/Make.rules
