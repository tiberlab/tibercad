# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))


all: tibercad

tibercad: $(objects)
	@$(CXX) $(CXXFLAGS) $(objects) -o $(topdir)/bin/tibercad \
		$(LIBS) $(LDFLAGS) -Wl,-export-dynamic -ltibercad  -lhetero -lifcore -limf -L /usr/pack/intel_fc-9.0-ma/lib/ -Wl,-rpath,/usr/pack/intel_fc-9.0-ma/lib
#	@$(CXX) $(CXXFLAGS) $(objects) $(TIBER_OBJECTS) 
#	-o $(topdir)/bin/tibercad $(LIBS) $(LDFLAGS)




include $(topdir)/Make.rules
