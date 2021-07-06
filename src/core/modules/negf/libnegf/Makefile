# $Id$

topdir := ../../../../..
include $(topdir)/Make.defs
include $(topdir)/Make.common

ARCH = $(shell uname -m)

CMAKE = cmake-3.18

CXXFLAGS = -fPIC -I . 


target := libnegf$(shared_libext)


libnegf$(shared_libext) :  .$(ARCH)/NegfWrapper$(obj-suffix) libnegf_$(ARCH).a
	@echo "Linking libnegf$(shared_libext) ..."
	$(CXX) -shared -o $@ $< \
		-Wl,--version-script,export.map \
		-Wl,-Bsymbolic -Wl,-Bstatic \
		-L. libnegf_$(ARCH).a \
		-Wl,-Bdynamic

.$(ARCH)/NegfWrapper$(obj-suffix) : libnegf_$(ARCH).a 

libnegf_$(ARCH).a: update
	@echo "Compiling libNEGF ...";
	cd libnegf; \
	MKLROOT=$(MKLROOT) FC=$(FC) $(CMAKE) . ; \
	$(MAKE) ; cd .. ; \
	cp libnegf/src/libnegf.a .	
	cp libnegf/src/api/fortran.h .
	cp libnegf/src/api/libnegf.hpp . 

update: checkout
	#cd libnegf; git pull origin; cd ..;

checkout:
	@echo "Checkout libNEGF..."
	@if [ ! -e libnegf ]; then git clone https://github.com/libnegf/libnegf.git; fi

clean:
	@echo "Clean libnegf..."
	rm -f libnegf_$(ARCH).a libnegf.so
	if test -d libnegf/src/api ; then make -C libnegf/src/api clean; fi	


include $(topdir)/Make.rules
