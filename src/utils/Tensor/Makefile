# $Id$

topdir := ../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.cpp)

#
# object files
objects		:= $(patsubst %.cpp, %.$(obj-suffix), $(srcfiles))



.PHONY: all
all:  $(topdir)/lib/libtensor.a


$(topdir)/lib/libtensor.a: $(objects)
	@ar cru $(topdir)/lib/libtensor.a $(objects)
	@ranlib $(topdir)/lib/libtensor.a

$(topdir)/lib/libtensor.$(shared_libext): $(objects)
	@$(CXX) -fPIC -shared -o $(topdir)/lib/libtensor.$(shared_libext) $(objects)




include $(topdir)/Make.rules
