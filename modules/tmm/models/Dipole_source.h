/*
 * incidentwave.h
 *
 *  Created on: 4 Oct 2021
 *      Author: pamiri
 */

#ifndef SRC_CORE_MODULES_TMM_MODELS_INCIDENTWAVE_H_
#define SRC_CORE_MODULES_TMM_MODELS_INCIDENTWAVE_H_

#include "TmmBoundaryModel.h"
#include "Tmm.h"

namespace libMesh
{
  class Elem;
}
class TBDLLOCAL Dipole_source : public TmmBoundaryModel
{

  public:

    //! Destructor
    ~Dipole_source(void) {};

    //! Creator function

    static Dipole_source* create(const ModelOptions& options);
    virtual void Calculate_M_Matrix(void);


  protected:

    virtual void do_init(void);




  private:
    double _kr;
    double _steps;
    //! Constructor
    Dipole_source(const ModelOptions& options);


};



inline
Dipole_source::Dipole_source(const ModelOptions& options) :
  TmmBoundaryModel(options)
{
}



inline
Dipole_source*
Dipole_source::create(const ModelOptions& options)
{
  return new Dipole_source(options);
}

#endif /* SRC_CORE_MODULES_TMM_MODELS_INCIDENTWAVE_H_ */
