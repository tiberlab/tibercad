# $Id: Makefile 415 2007-07-17 09:35:34Z maufder $

topdir := ../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.cpp)

#
# object files
objects		:= $(patsubst %.cpp, %.o, $(srcfiles))



.PHONY: all clean


all: $(objects)
	@$(CXX) -fPIC -shared -o $(topdir)/lib/libnewmat.so $(objects)

#@ar rcs $(topdir)/lib/libnewmat.a $(objects)



%.o: %.cpp
	@$(CXX) -fPIC $(TIBER_INCLUDES) -c $<


include $(topdir)/Make.rules
