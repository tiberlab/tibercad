# $Id$

topdir := ../../..

include $(topdir)/Make.common

target := .$(ARCH)/liblicense.a

ifeq ($(enable_license_check), yes)
  CPPFLAGS += -I/home/maufder/Projects/TIBERCAD/LM-X/lmx-sdk-4.9.8/include
  LMXCLIENTLIB := /home/maufder/Projects/TIBERCAD/LM-X/lmx-sdk-4.9.8/linux_x64/liblmxclient.a

.$(ARCH)/liblicense.a: $(objects)
	@echo "Creating $@ ..."
	@cp $(LMXCLIENTLIB) $@
	@$(AR) ru $@ $(objects)
	@$(RANLIB) $@
else
.$(ARCH)/liblicense.a: $(objects)
	@echo "Creating $@ ..."
	$(AR) cr $@ $(objects)
	@$(RANLIB) $@
endif



distclean:
	@rm -f .$(ARCH)/liblicense.a 

include $(topdir)/Make.rules
