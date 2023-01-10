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
class TBDLLOCAL Incident_wave : public TmmBoundaryModel
{

  public:

    //! Destructor
    ~Incident_wave(void) {};

    //! Creator function


    static Incident_wave* create(const ModelOptions& options);
    virtual void Calculate_M_Matrix(void);


  protected:

    virtual void do_init(void);




  private:

    //! Constructor
    Incident_wave(const ModelOptions& options);


};



inline
Incident_wave::Incident_wave(const ModelOptions& options) :
  TmmBoundaryModel(options)
{
}



inline
Incident_wave*
Incident_wave::create(const ModelOptions& options)
{
  return new Incident_wave(options);
}

#endif /* SRC_CORE_MODULES_TMM_MODELS_INCIDENTWAVE_H_ */
