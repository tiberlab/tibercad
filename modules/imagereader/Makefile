# $Id: Makefile 3542 2013-03-01 09:31:59Z maufder $

topdir := ../../../..

include $(topdir)/Make.common

target := modules
#target :=

modules: module.conf paintlib-2.6.2
	@$(BUILDMODULE) module.conf
	@cp -L paintlib-2.6.2/lib/.libs/libpaintlib.so.1 $(topdir)/$(ARCH)/modules/imagereader


paintlib-2.6.2:
	@tar zxf paintlib-2.6.2.tar.gz
	@cd paintlib-2.6.2 ; \
	 patch -p0 < ../paintlib.patch; \
	 CC=$(CC) CXX=$(CXX) CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" \
	 ./configure && make


include $(topdir)/Make.rules
