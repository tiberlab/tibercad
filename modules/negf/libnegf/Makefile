# $Id$

topdir := ../../..
include $(topdir)/Make.defs
include $(topdir)/Make.common

ARCH = $(shell uname -m)

CMAKE = cmake -DMKL_LIBDIR=/opt/intel/oneapi/mkl/latest/lib/intel64 -DLAPACK_LIBRARY_DIR=/opt/intel/oneapi/mkl/latest/lib/intel64\;/usr/lib/x86_64-linux-gnu -DSCALAPACK_LIBRARY_DIR=/opt/intel/oneapi/mkl/latest/lib/intel64

CXXFLAGS = -fPIC -I . 

CUDADIR = /usr/pack/cudatoolkit-11.4-ma/lib64/
CUDALIBS = -L$(CUDADIR) -Wl,-rpath=$(CUDADIR) -lcudart -lcublas -lcusolver -lstdc++

target := libnegf$(shared_libext)


libnegf$(shared_libext) :  .$(ARCH)/NegfWrapper$(obj-suffix) libnegf_$(ARCH).a
	@echo "Linking libnegf$(shared_libext) ..."
	dep=$$(find . -name libmpifx.a); \
	mpi_link=$$(mpifort -showme:link 2>/dev/null || true); \
	$(CXX) -shared -o $@ $< \
		-Wl,--version-script,export.map \
		-Wl,-Bsymbolic -Wl,-Bstatic \
		-L. libnegf_$(ARCH).a $$dep \
		-Wl,-Bdynamic \
		$$mpi_link -lmpi_usempif08 -lmpi -lgfortran

.$(ARCH)/NegfWrapper$(obj-suffix) : libnegf_$(ARCH).a 


#cp ../myCUDACompiler.cmake CMakeFiles/3.18.4/CMakeCUDACompiler.cmake;
libnegf_$(ARCH).a: update
	@echo "Compiling libNEGF ...";
	mkdir -p build; cd build; \
	$(CMAKE) -C ../config.cmake -S ../libnegf -B .; \
	$(MAKE); cd .. ;
	cp build/src/libnegf.a libnegf_$(ARCH).a	
	cp build/src/api/libnegf.hpp . 
	cp libnegf/src/api/fortran.h .
	cp libnegf/src/api/lnParams.h .

update: checkout
	cd libnegf; git pull origin; cd ..;

checkout:
	@echo "Checkout libNEGF..."
	@if [ ! -e libnegf ]; then \
		git clone https://github.com/libnegf/libnegf.git; \
		cd libnegf; \
	fi

clean:
	@echo "Clean libnegf..."
	rm -f libnegf_$(ARCH).a libnegf.so
	if test -d build/src/api ; then make -C build/src/api clean; fi	


include $(topdir)/Make.rules
