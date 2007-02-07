# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))



all: $(objects)


tibercad: all
	@$(CXX) $(CXXFLAGS) $(objects) -o $(topdir)/bin/tibercad \
		$(LIBS) $(LDFLAGS) -ltibercad
#	@$(CXX) $(CXXFLAGS) $(objects) $(TIBER_OBJECTS) \
#	-o $(topdir)/bin/tibercad $(LIBS) $(LDFLAGS)




include $(topdir)/Make.rules
