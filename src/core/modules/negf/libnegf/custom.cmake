# Selecting toolchain specific settings via environment variable, so that multiple shells with
#  different toolchains are possible in parallel.
set(TOOLCHAIN "gnu" CACHE STRING
  "Selects architecture and architecture dependant settings to be used")

# Get host name, so that settings can be host dependent
site_name(host_name)

set(FC "gfortran-11.1" CACHE STRING "Fortran compiler")
set(CC "gcc-11.1" CACHE STRING "C compiler")
set(MYINCLUDE "$ENV{HOME}/include/gnu-11.1.0" CACHE STRING "C compiler")

#
# Overriding (pre-populating) options in config.cmake
#
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Release|Debug)")

option(BUILD_API "Whether API should be built" FALSE)

option(WITH_MPI "Whether MPI should be used" TRUE)

option(WITH_OMP "Whether OMP should be used" TRUE)

option(WITH_SOCKETS "Whether socket communication should allowed for" FALSE)

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)

option(WITH_TRANSPORT "Whether transport via libNEGF should be enabled" TRUE)
option(WITH_TRANSPORT_GPU "Whether libNEGF should enable GPU" TRUE)

if(NOT WITH_MPI)
  option(WITH_ARPACK "Whether the ARPACK library should be included" FALSE)
else()
  option(WITH_ARPACK "Whether the ARPACK library should be included" FALSE)
endif()

if(WITH_MPI)
  option(WITH_ELSI "Whether the ELSI libraries are compiled" FALSE)
  option(ELSI_WITH_PEXSI "Whether the ELSI libraries are compiled  with PEXSI support" FALSE)
endif()

set(MKL_LIBDIR "/usr/pack/intel_mkl-11.1-ma/mkl/lib/intel64" CACHE STRING
	"Directory where to look for mkl libraries"	)
set(INTEL_DIR "/usr/pack/intel_mkl-11.1-ma//lib/intel64" CACHE STRING
	"Directory where to look for intel libraries"	)
#set(GNU_DIR "/usr/pack/gcc-7.2.0-ma/lib64" CACHE STRING
#	"Directory where to look for gnu libraries"	)
set(IFORT_DIR "/usr/pack/intel_ifort-13.0-ma/lib/intel64" CACHE STRING
	"Directory where to look for ifort libraries"	)
set(ENV{MKLROOT} ${MKL_LIBDIR})
#set(ENV{CMAKE_PREFIX_PATH} ${MKL_LIBDIR})

if (WITH_MPI)
   set(MPI_DIR "/usr/pack/mpich-3.1.1-ma/lib" CACHE STRING
	   "Directory where to look for MPI libs"	)
endif()
#
# Overriding (pre-populating) options of the toolchain file
#

if("${TOOLCHAIN}" STREQUAL "gnu")

  set(TEST_OMP_THREADS "2" CACHE STRING "Number of threads to use for each test")
  set(TEST_MPI_PROCS "2" CACHE STRING "Number of mpi processes to use for each test")


  if(WITH_MPI)

    if(WITH_TRANSPORT_GPU)
      set(CMAKE_CUDA_COMPILER "/usr/pack/cudatoolkit-11.4-ma/bin/nvcc" CACHE STRING "CUDA-C compiler")
    endif()

    set(CMAKE_Fortran_COMPILER "mpif90"
       CACHE STRING "Fortran compiler")
  
    set(CMAKE_Fortran_FLAGS "-fc=${FC} -fmax-errors=5 -fPIC -I${MYINCLUDE}" 
       CACHE STRING "General Fortran flags")

    set(CMAKE_C_COMPILER "mpicc" CACHE STRING "C compiler")

    set(CMAKE_C_FLAGS "-cc=${CC} -fmax-errors=5 -fPIC -I${MYINCLUDE}" 
       CACHE STRING "General Fortran flags")

    set(SCALAPACK_LIBRARY_DIR "${MKL_LIBDIR}"
    	   CACHE STRING "Directories where Scalapack libraries can be found")

    set(SCALAPACK_LIBRARY "mkl_scalapack_lp64;mkl_blacs_intelmpi_lp64" CACHE STRING 
    	   "Scalapack libraries to link")

    set(TEST_RUNNER_TEMPLATE
    	   "env OMP_NUM_THREADS=\${TEST_OMP_THREADS} mpiexec --bind-to none -n \${TEST_MPI_PROCS}"
    	   CACHE STRING "How to run the tests")

    set(OTHER_LIBRARY_DIRS ""
    		CACHE STRING "Directories where the other libraries can be found")
    set(OTHER_INCLUDE_DIRS "" CACHE STRING "Other include directories to consider")

  else()
    set(CMAKE_Fortran_COMPILER "gfortran-7.2.0" CACHE STRING "Fortran compiler")
    set(CMAKE_C_COMPILER "gcc-7.2.0" CACHE STRING "C compiler")
  endif()

  set(LAPACK_LIBRARY_DIR "${MKL_LIBDIR};/usr/lib/x86_64-linux-gnu"
      CACHE STRING "Paths to scan when looking for the LAPACK/BLAS libraries")
      
  set(LAPACK_LIBRARY "mkl_intel_lp64;mkl_gnu_thread;mkl_core;-lgomp" CACHE STRING 
    	    "LAPACK and BLAS libraries to link")
     
        #set(ELSI_ROOT "$ENV{HOME}/opt/elsi/gnu8-omp3.1" CACHE STRING
        #  "Root directory of the ELSI installation")
      
        #set(PEXSI_EXTERNAL_LIBRARY_DIRS "$ENV{HOME}/opt/gcc/8.3/lib64" CACHE STRING
        #  "Directories with PEXSI libraries")
      
else()
      
  message(FATAL_ERROR "Uknown toolchain '$ENV{TOOLCHAIN}'")

endif()

