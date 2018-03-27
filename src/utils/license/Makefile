# $Id$

topdir := ../../..

include $(topdir)/Make.common

target := .$(ARCH)/liblicense.a

ifeq ($(enable_license_check), yes)
  CPPFLAGS += -I/home/maufder/Projects/TIBERCAD/LM-X/lmx-sdk/include
  LMXCLIENTLIB := /home/maufder/Projects/TIBERCAD/LM-X/lmx-sdk/linux_x64/liblmxclient.a

.$(ARCH)/liblicense.a: $(objects)
	@echo "Creating $@ ..."
	@cp $(LMXCLIENTLIB) $@
	@$(AR) ru $@ $(objects)
	@$(RANLIB) $@
else
.$(ARCH)/liblicense.a: $(objects)
	@echo "Creating $@ ..."
	@$(AR) cu $@ $(objects)
	@$(RANLIB) $@
endif



distclean:
	@rm -f .$(ARCH)/liblicense.a 

include $(topdir)/Make.rules
