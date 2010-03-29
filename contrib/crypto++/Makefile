
topdir := ../..

include $(topdir)/Make.common

#package := cryptopp552.zip
package := cryptopp560.zip

all: build/libcryptopp.a

build/libcryptopp.a: stamp-h
	@$(MAKE) CC=$(CC) CXX=$(CXX) -C build


stamp-h: $(package)
	@rm -rf build ; mkdir build && unzip $(package) -d build/ && touch stamp-h

clean:
	@rm -rf build ; rm stamp-h

distclean: clean

include $(topdir)/Make.rules
