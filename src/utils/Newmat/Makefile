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


all: $(objects)
	@$(CXX) -fPIC -shared -o $(topdir)/lib/libnewmat.so $(objects)





include $(topdir)/Make.rules
