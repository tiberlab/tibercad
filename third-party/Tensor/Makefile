# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.cpp)

#
# object files
objects		:= $(patsubst %.cpp, .$(ARCH)/%$(obj-suffix), $(srcfiles))

target := $(ARCH)/libtensor.a



$(ARCH)/libtensor.a: $(objects)
	@mkdir -p $(ARCH)
	@rm -f $@
	@$(AR) cr $@ $(objects)
	@$(RANLIB) $@

$(libdir)/libtensor$(shared_libext): $(objects)
	@$(CXX) -fPIC -shared -o $@ $(objects)




include $(topdir)/Make.rules
