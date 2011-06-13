# $Id$

topdir := ../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.cpp)

#
# object files
objects		:= $(patsubst %.cpp, .$(ARCH)/%$(obj-suffix), $(srcfiles))

target := $(libdir)/libtensor.a



$(libdir)/libtensor.a: $(objects)
	@$(AR) cru $@ $(objects)
	@$(RANLIB) $@

$(libdir)/libtensor$(shared_libext): $(objects)
	@$(CXX) -fPIC -shared -o $@ $(objects)




include $(topdir)/Make.rules
