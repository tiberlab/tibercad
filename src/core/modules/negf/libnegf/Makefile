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

.$(ARCH)/NegfWrapper$(obj-suffix) : libnegf_$(ARCH).a 

libnegf_$(ARCH).a: update
	@echo "Updating libNEGF ..."  
	cd libnegf; git status; \
	echo "Compiling libNEGF ...";\
	$(MAKE) -C src/api ; cd .. ; 
	mv libnegf/src/api/libnegf_$(ARCH).a .	
	cp libnegf/src/api/fortran.h .
	cp libnegf/src/api/libnegf.h . 

update: checkout
	#cd libnegf; git pull origin; cd ..;

checkout:
	@echo "Checkout libNEGF..."
	@if [ ! -e libnegf ]; then git clone https://pecchia@bitbucket.org/pecchia/libnegf.git; fi

clean:
	@echo "Clean libnegf..."
	rm -f libnegf_$(ARCH).a libnegf.so
	cd libnegf/src/api; make clean;	


include $(topdir)/Make.rules
