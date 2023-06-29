# $Id: Makefile.in 2839 2011-06-26 14:52:41Z maufder $

topdir := ../..

include $(topdir)/Make.common

ifeq ($(target_os), linux)
 PACKAGEFLAGS := -Wl,-rpath,'$$ORIGIN/../$(relative_libdir)' -Wl,-rpath-link,$(libdir)
endif
ifeq ($(target_os), darwin)
 PACKAGEFLAGS := -Wl,-rpath,@loader_path/.. -Wl,-rpath,@loader_path/../extern 
endif



srcfiles := build_module.C
objects	:= $(patsubst %.C, .$(ARCH)/%$(obj-suffix), $(srcfiles))
objects += $(topdir)/src/base/common/.$(ARCH)/Utils.$(METHOD).o
objects += $(topdir)/src/base/common/.$(ARCH)/ModelOptions.$(METHOD).o
objects += $(topdir)/src/base/io/.$(ARCH)/InputParser.$(METHOD).o
objects += $(topdir)/src/base/io/.$(ARCH)/Messages.$(METHOD).o

cross_srcfiles := $(srcfiles)
cross_srcfiles += $(topdir)/src/base/common/Utils.C
cross_srcfiles += $(topdir)/src/base/common/ModelOptions.C
cross_srcfiles += $(topdir)/src/base/io/InputParser.C
cross_srcfiles += $(topdir)/src/base/io/Messages.C
BUILD_CPPFLAGS := -DBUILD=\"$(BUILD)\" $(CPPFLAGS)
BUILD_CXXFLAGS := $(CXXFLAGS)
BUILD_LDFLAGS := -Wl,-rpath,'$$ORIGIN/../extern/$(BUILD)/lib' -L../../extern/$(BUILD)/lib -lboost_regex -lboost_filesystem -lboost_system
BUILD_CXX := g++


SUFF := 
ifeq ($(target_os),win32)
  SUFF := .exe
  #objects += tibercad_res.o
  WINDRES := $(subst gcc,windres,$(CC))
  LDFLAGS += -lcomdlg32
endif

ifeq ($(target_os),darwin)
  ASNEEDED := 
else
  ASNEEDED := -Wl,--as-needed
endif


all: $(topdir)/$(ARCH)/bin/build_module$(SUFF)

ifeq ($(cross_compiling), yes)
all: $(topdir)/bin/build_module-$(BUILD)
endif

$(topdir)/$(ARCH)/bin/build_module$(SUFF): $(objects) 
	@echo "Compiling build_module executable ..."
	@mkdir -p $(topdir)/$(ARCH)/bin
	@$(CXX) $(CXXFLAGS) $(objects) $(ASNEEDED) $(LDFLAGS) -o $(topdir)/$(ARCH)/bin/build_module$(SUFF) $(BOOST_LIBS)
	@mkdir -p $(topdir)/SDK/$(ARCH)/bin
	@cp $(topdir)/$(ARCH)/bin/build_module$(SUFF) $(topdir)/SDK/$(ARCH)/bin
	
	
#	@$(CXX) $(CXXFLAGS) $(objects) -Wl,--as-needed $(LDFLAGS) -o $(topdir)/$(ARCH)/bin/build_module$(SUFF) $(BOOST_LIBS)
#	@$(CXX) $(CXXFLAGS) $(objects) -o $(topdir)/$(ARCH)/bin/build_module$(SUFF) $(BOOST_LIBS) $(PACKAGEFLAGS) $(LDFLAGS)

$(topdir)/bin/build_module-$(BUILD): $(cross_srcfiles)
	@echo "compiling build_module for cross compilation ..."
	@$(BUILD_CXX) $(BUILD_CPPFLAGS) $(BUILD_CXXFLAGS) $(BUILD_LDFLAGS) \
		-o $@ $(BUILD_LDFLAGS) $(cross_srcfiles)




tibercad_res.o : $(topdir)/share/tibercad.rc $(topdir)/share/tibercad.ico $(topdir)/include/base/common/tiber_config.h
	$(WINDRES) -I $(topdir)/share -I $(topdir)/include/base/common \
	  -O coff -D SVNREVISION=$(SVNREVISION) $< -o $@


include $(topdir)/Make.rules
