# $Id$

topdir := ../../..

include $(topdir)/Make.common

CPPFLAGS += -I$(topdir)/contrib/crypto++

#
# source files
srcfiles 	:= License.$(obj-suffix)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))



.PHONY: all
all:  liblicense.a create_license

# we take the library of crypto++ and add our stuff...
liblicense.a: $(topdir)/contrib/crypto++/libcryptopp.a $(objects)
	@cp $(topdir)/contrib/crypto++/libcryptopp.a $@
	@ar ru $@ $(objects)
	@ranlib $@


create_license: create_license.$(obj-suffix) liblicense.a
	$(CXX) -o $@ $< liblicense.a -lpthread


$(topdir)/contrib/crypto++/libcryptopp.a: $(topdir)/contrib/crypto++/stamp-h
	@$(MAKE) -C $(topdir)/contrib/crypto++


$(topdir)/contrib/crypto++/stamp-h: $(topdir)/contrib/crypto++/cryptopp552.zip
	@cd $(topdir)/contrib/crypto++ && unzip cryptopp552.zip && touch stamp-h


distclean:
	@rm -f liblicense.a create_license

include $(topdir)/Make.rules
