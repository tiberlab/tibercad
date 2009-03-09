// $Id$

#ifndef _EIGENVALUEPROBLEM_H_
#define _EIGENVALUEPROBLEM_H_

#include "SimulationInterface.h"


//! Abstract class to solve complex valued eigenvalue problem
class EigenvalueProblem : public SimulationInterface
{

  public:

    //! Constructor
    EigenvalueProblem(void) { };

    //! Destructor
    ~EigenvalueProblem(void) { };


  protected:


  private:

};


#endif // _EIGENVALUEPROBLEM_H_
