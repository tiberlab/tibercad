// $Id: PermittivityModel.h 4729 2018-12-05 07:58:16Z maufder $

#ifndef _INCOHERENCE_H_
#define _INCOHERENCE_H_

#include "Tmm.h"
#include "TmmBulkModel.h"

namespace libMesh
{
  class Elem;
}

// Base class for InCoherence model
class TBDLEXPORT IncoherentModel : public TmmBulkModel
{

  public:

  virtual ~IncoherentModel(void) {};
  
  const double& get_Incoherent_Index(void) const;
  static IncoherentModel* create(const ModelOptions& options);

  
protected:
  
    IncoherentModel(const ModelOptions& options);

    void set_Incoherent_Index(const double& Incoheret_Index);




  private:

   double _Incoheret_Index;

};




#endif // _POLARIZATIONMODEL_H_
