CFLAGS	   = 
FFLAGS     =
CPPFLAGS   = 
FPPFLAGS   =
LOCDIR     = ./

SLEPC_DIR  = /usr/pack/slepc-2.2.1-ma
PETSC_DIR  = /usr/pack/petsc-2.2.1-ma

BOPT = O_complex

MANSEC           = EPS

include ${SLEPC_DIR}/bmake/slepc_common


all: ex7.o 
	-${CLINKER}  -o ex7 $(CXXFLAGS) $(FFLAGS) $(INCLUDE) ex7.o ${LDFLAGS} ${SLEPC_LIB}
	\cp ex7 ../../../bin/





