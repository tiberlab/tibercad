# $Id$

topdir := ../../../../..
include $(topdir)/Make.defs
include $(topdir)/Make.common

ARCH = $(shell uname -m)

CMAKE = cmake

CXXFLAGS = -fPIC -I . 


target := libnegf$(shared_libext)


libnegf$(shared_libext) :  .$(ARCH)/NegfWrapper$(obj-suffix) libnegf_$(ARCH).a
	@echo "Linking libnegf$(shared_libext) ..."
	dep=$$(find . -name libmpifx.a); \
	$(CXX) -shared -o $@ $< \
		-Wl,--version-script,export.map \
		-Wl,-Bsymbolic -Wl,-Bstatic \
		-L. libnegf_$(ARCH).a $$dep \
		-Wl,-Bdynamic

.$(ARCH)/NegfWrapper$(obj-suffix) : libnegf_$(ARCH).a 

libnegf_$(ARCH).a: update
	@echo "Compiling libNEGF ...";
	mkdir -p build; cd build; \
	$(CMAKE) -C ../config.cmake -S ../libnegf -B .; \
	$(MAKE) ; cd .. ;
	cp build/src/libnegf.a libnegf_$(ARCH).a	
	cp build/src/api/libnegf.hpp . 
	cp libnegf/src/api/fortran.h .
	cp libnegf/src/api/lnParams.h .

update: checkout
	#cd libnegf; git pull origin; cd ..;

checkout:
	@echo "Checkout libNEGF..."
	@if [ ! -e libnegf ]; then \
		git clone https://github.com/libnegf/libnegf.git; \
		cd libnegf; \
		patch -p1 < ../libnegf.patch; \
	fi

clean:
	@echo "Clean libnegf..."
	rm -f libnegf_$(ARCH).a libnegf.so
	if test -d build/src/api ; then make -C build/src/api clean; fi	


include $(topdir)/Make.rules
