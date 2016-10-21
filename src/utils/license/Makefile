# $Id$

topdir := ../../..

include $(topdir)/Make.common

CPPFLAGS += -I$(topdir)/contrib/crypto++/$(ARCH)

target := .$(ARCH)/liblicense.a $(topdir)/bin/create_license


# we take the library of crypto++ and add our stuff...
.$(ARCH)/liblicense.a: $(topdir)/contrib/crypto++/$(ARCH)/libcryptopp.a $(objects)
	@echo "Creating $@ ..."
	@cp $(topdir)/contrib/crypto++/$(ARCH)/libcryptopp.a $@
	@$(AR) ru $@ $(objects)
	@$(RANLIB) $@


$(topdir)/bin/create_license: .$(ARCH)/create_license$(obj-suffix) .$(ARCH)/liblicense.a
ifeq ($(target_os),win32)
	$(CXX) -o $@ $< .$(ARCH)/liblicense.a
else
	$(CXX) -o $@ $< .$(ARCH)/liblicense.a $(LDFLAGS) -lpthread
endif


$(topdir)/contrib/crypto++/$(ARCH)/libcryptopp.a:
	@$(MAKE) -C $(topdir)/contrib/crypto++



distclean:
	@rm -f .$(ARCH)/liblicense.a $(topdir)/bin/create_license
	@$(MAKE) -C $(topdir)/contrib/crypto++ clean


include $(topdir)/Make.rules
