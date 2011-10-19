# $Id: Makefile.in 2839 2011-06-26 14:52:41Z maufder $

topdir := ../..

include $(topdir)/Make.common


srcfiles := build_module.C
objects	:= $(patsubst %.C, .$(ARCH)/%$(obj-suffix), $(srcfiles))
objects += $(topdir)/src/base/common/.$(ARCH)/Utils.$(METHOD).o
objects += $(topdir)/src/base/common/.$(ARCH)/ModelOptions.$(METHOD).o
objects += $(topdir)/src/base/io/.$(ARCH)/InputParser.$(METHOD).o
objects += $(topdir)/src/base/io/.$(ARCH)/Messages.$(METHOD).o



SUFF := 
ifeq ($(target_os),win32)
  SUFF := ".exe"
  #objects += tibercad_res.o
  WINDRES := $(subst gcc,windres,$(CC))
  #ADDITIONAL_LIBS += -lcomdlg32
endif

ifeq ($(target_os),darwin)
  ASNEEDED := 
else
  ASNEEDED := "-Wl,--as-needed"
endif



all: build_module

build_module: $(objects) 
	@echo "Compiling build_module executable ..."
	@$(CXX) $(CXXFLAGS) $(objects) -o $(topdir)/$(ARCH)/bin/build_module$(SUFF) $(ASNEEDED) $(LDFLAGS)





tibercad_res.o : $(topdir)/share/tibercad.rc $(topdir)/share/tibercad.ico $(topdir)/include/base/common/tiber_config.h
	$(WINDRES) -I $(topdir)/share -I $(topdir)/include/base/common \
	  -O coff -D SVNREVISION=$(SVNREVISION) $< -o $@


include $(topdir)/Make.rules
