# $Id$

topdir := ../../../..

MODULE_NAME := excitontransport

include $(topdir)/Make.module


MODULE_LIB =  $(module_dir)/$(MODULE_NAME)/excitontransport$(shared_libext)

$(eval $(call MODULE, ExcitonTransport.C ExcitonProperties.C))

$(eval $(call MODULE, ExcitonModel.C, $(MODULE_LIB)))


include $(topdir)/Make.rules
