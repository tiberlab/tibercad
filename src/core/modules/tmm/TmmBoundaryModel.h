// $Id$

#ifndef _TMMBOUNDARYMODEL_H_
#define _TMMBOUNDARYMODEL_H_

#include "PhysicalModel.h"


namespace libMesh
{
  class Elem;
  class Point;
}

//! The base class for thermal balance boundary conditions
class TmmBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    ~TmmBoundaryModel(void) {};

    //! Creator function
    static TmmBoundaryModel* create(const MaterialBoundary* boundary,
        const ModelOptions& options);





  protected:

    //! Constructor
    TmmBoundaryModel(const ModelOptions& options);



  private:


};


inline
TmmBoundaryModel::TmmBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}



#endif // _THERMALBOUNDARYMODEL_H_
