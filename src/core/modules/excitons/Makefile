# $Id$

topdir := ../../../..

MODULE_NAME := exciton

include $(topdir)/Make.module

MODULELIB :=  $(topdir)/lib/modules/$(MODULE_NAME)/excitontransport$(shared_libext)



$(eval $(call MODULE, ExcitonTransport.C ExcitonProperties.C))
$(eval $(call MODULE, ExcitonModel.C, $(MODULELIB)))


include $(topdir)/Make.rules
