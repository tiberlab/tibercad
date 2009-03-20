# $Id$

topdir := ../../..

include $(topdir)/Make.common

CPPFLAGS += -I$(topdir)/contrib/crypto++

target := liblicense.a $(topdir)/bin/create_license


# we take the library of crypto++ and add our stuff...
liblicense.a: $(topdir)/contrib/crypto++/libcryptopp.a $(objects)
	@cp $(topdir)/contrib/crypto++/libcryptopp.a $@
	@ar ru $@ $(objects)
	@ranlib $@


$(topdir)/bin/create_license: create_license$(obj-suffix) liblicense.a
	$(CXX) -o $@ $< liblicense.a -lpthread


$(topdir)/contrib/crypto++/libcryptopp.a: $(topdir)/contrib/crypto++/stamp-h
	@$(MAKE) CC=$(CC) CXX=$(CXX) -C $(topdir)/contrib/crypto++


$(topdir)/contrib/crypto++/stamp-h: $(topdir)/contrib/crypto++/cryptopp552.zip
	@cd $(topdir)/contrib/crypto++ && unzip cryptopp552.zip && touch stamp-h


distclean:
	@rm -f liblicense.a $(topdir)/bin/create_license
	@$(MAKE) -C $(topdir)/contrib/crypto++ clean

include $(topdir)/Make.rules
