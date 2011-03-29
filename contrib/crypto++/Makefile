
topdir := ../..

include $(topdir)/Make.common

#package := cryptopp552.zip
package := cryptopp560.zip

all: $(ARCH)/libcryptopp.a

$(ARCH)/libcryptopp.a: stamp-h
	@$(MAKE) CC=$(CC) CXX=$(CXX) -C $(ARCH)


stamp-h: $(package)
	@rm -rf $(ARCH) ; mkdir $(ARCH) && unzip $(package) -d $(ARCH)/ && touch stamp-h

clean:
	@rm -rf $(ARCH)
	@-rm stamp-h

distclean: clean

include $(topdir)/Make.rules
