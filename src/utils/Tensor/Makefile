# $Id$

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
	@$(CXX) -fPIC -shared -o $(topdir)/lib/libtensor.so $(objects)

#@ar rcs $(topdir)/lib/libtensor.a $(objects)



%.o: %.cpp
	@$(CXX) -fPIC $(TIBER_INCLUDES) -c $<


include $(topdir)/Make.rules
