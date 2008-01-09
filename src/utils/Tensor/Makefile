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


all: $(objects)
	@$(CXX) -fPIC -shared -o $(topdir)/lib/libtensor.$(shared_libext) $(objects)

#@ar rcs $(topdir)/lib/libtensor.a $(objects)



include $(topdir)/Make.rules
