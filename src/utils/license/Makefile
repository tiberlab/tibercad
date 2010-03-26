# $Id$

topdir := ../../..

include $(topdir)/Make.common

CPPFLAGS += -I$(topdir)/contrib/crypto++/build

target := liblicense.a $(topdir)/bin/create_license


# we take the library of crypto++ and add our stuff...
liblicense.a: $(topdir)/contrib/crypto++/build/libcryptopp.a $(objects)
	@cp $(topdir)/contrib/crypto++/build/libcryptopp.a $@
	@ar ru $@ $(objects)
	@ranlib $@


$(topdir)/bin/create_license: create_license$(obj-suffix) liblicense.a
	$(CXX) -o $@ $< liblicense.a -lpthread


$(topdir)/contrib/crypto++/build/libcryptopp.a:
	@$(MAKE) -C $(topdir)/contrib/crypto++



distclean:
	@rm -f liblicense.a $(topdir)/bin/create_license
	@$(MAKE) -C $(topdir)/contrib/crypto++ clean


include $(topdir)/Make.rules
