# $Id$

topdir := ../../../../..
include $(topdir)/Make.defs
ARCH = $(shell uname -m)

CXXFLAGS = -fPIC -I . 


target := libnegf$(shared_libext)


libnegf$(shared_libext) :  .$(ARCH)/NegfWrapper$(obj-suffix) libnegf_$(ARCH).a
	@echo "Linking libnegf$(shared_libext) ..."
	$(CXX) -shared -o $@ $< \
		-Wl,--version-script,export.map \
		-Wl,-Bsymbolic -Wl,-Bstatic \
		-L. libnegf_$(ARCH).a \
		-Wl,-Bdynamic



include $(topdir)/Make.rules
