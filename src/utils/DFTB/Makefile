CFLAGS	   = 
CXXLAGS	   = 
FFLAGS     =
CPPFLAGS   = 
FPPFLAGS   =

DFTB_DIR  =.

CXX        = g++-4.1.1


CXXFLAGS += -fPIC $(CC_INCLUDES) -I ../../../include/utils/DFTB


all: libdftbpc.so


libdftbpc.so:  DftbpWrapper.o
	@$(CXX) -shared -Wl,-soname,libdftbpc.so \
		-Wl,--version-script,export.map \
		-o libdftbpc.so $<  \
		-Wl,-Bsymbolic -Wl,-Bstatic \
		-L. libdftbp.a \
		-Wl,-Bdynamic
	@cp libdftbpc.so ../../../lib

%.o: %.C
	@$(CXX) $(CXXFLAGS) -c $<


clean:
	rm -f *.o libdftbpc.so

distclean: clean



#include $(SLEPC_DIR)/bmake/slepc_common_variables
