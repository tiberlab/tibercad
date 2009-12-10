// $Id$


#ifndef _TIBERPETSCUTILS_H_
#define _TIBERPETSCUTILS_H_

#include "PetscRuntimeError.h"

#include "petscksp.h"


class ModelOptions;

//! Provides utilities for PETSc
class TiberPetscUtils
{

  public:

    //! Extract KSPType from options
    static KSPType extract_KSPType(const ModelOptions& options);


    //! Extract PCType from options
    static PCType extract_PCType(const ModelOptions& options);


    //! Extract line search type (none, cubic, quadratic) from options
    static int extract_LSType(const ModelOptions& options);


    //! Check the PETSc error number
    /*!
     * Throws an exception if errorcode contains a PETSc error number
     */
    static void _checkerr(int errorcode, int line, const char* file);
#ifdef DEBUG
#define checkerr(errcode) _checkerr(errcode, __LINE__, __FILE__);
#else
     static void checkerr(int errorcode)
     {
       if (errorcode != 0)
         throw PetscRuntimeError(errorcode);
     }
#endif


  private:

    TiberPetscUtils();

};


//
// inline members
//



#endif // _TIBERPETSCUTILS_H_
