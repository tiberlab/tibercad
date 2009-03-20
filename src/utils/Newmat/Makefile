# $Id: Makefile 415 2007-07-17 09:35:34Z maufder $

topdir := ../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.cpp)

#
# object files
objects		:= $(patsubst %.cpp, %$(obj-suffix), $(srcfiles))


target := $(libdir)/libnewmat.a


$(libdir)/libnewmat.a: $(objects)
	@$(AR) cru $@ $(objects)
	@ranlib $@

$(libdir)/libnewmat$(shared_libext): $(objects)
	@$(CXX) -fPIC -shared -o $@ $(objects)




include $(topdir)/Make.rules
