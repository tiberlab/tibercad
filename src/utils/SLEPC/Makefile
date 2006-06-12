CFLAGS	   = 
FFLAGS     =
CPPFLAGS   = 
FPPFLAGS   =
LOCDIR     =./
PETSC_ARCH = linux-gnu-complex-opt

SLEPC_DIR  = /usr/pack/slepc-2.3.0-ma
PETSC_DIR  = /usr/pack/petsc-2.3.0-ma

BOPT = O_complex

MANSEC           = EPS

include ${SLEPC_DIR}/bmake/slepc_common


all: ex7.o 
	@${CLINKER}  -o ex7 $(CXXFLAGS) $(FFLAGS) $(INCLUDE) ex7.o ${LDFLAGS} ${SLEPC_LIB}
	\cp ex7 ../../../bin/


clean:
	rm -f *.o ex7



