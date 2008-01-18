# $Id: Makefile 415 2007-07-17 09:35:34Z maufder $

topdir := ../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.cpp)

#
# object files
objects		:= $(patsubst %.cpp, %.$(obj-suffix), $(srcfiles))



.PHONY: all
all: $(topdir)/lib/libnewmat.a


$(topdir)/lib/libnewmat.a: $(objects)
	@$(AR) cru $(topdir)/lib/libnewmat.a $(objects)
	@ranlib $(topdir)/lib/libnewmat.a

$(topdir)/lib/libnewmat.$(shared_libext): $(objects)
	@$(CXX) -fPIC -shared -o $(topdir)/lib/libnewmat.$(shared_libext) $(objects)




include $(topdir)/Make.rules
