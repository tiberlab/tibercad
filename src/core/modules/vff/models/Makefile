# $Id: Makefile 2039 2010-07-22 16:00:20Z gromano $

topdir := ../../../../..


# we have to find the necessary include files of the module
MODULE_INCLUDES = -I..

include $(topdir)/Make.module

MODULE_LIB := $(module_dir)/$(MODULE_NAME)/vff$(shared_libext)
$(foreach lib, $(srcfiles), $(eval $(call MODULE, $(lib), $(MODULE_LIB))))


include $(topdir)/Make.rules
