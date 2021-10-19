/*
 * Mirror.h
 *
 *  Created on: 30 Sep 2021
 *      Author: pamiri
 */

#ifndef SRC_CORE_MODULES_TMM_MODELS_MIRROR_H_
#define SRC_CORE_MODULES_TMM_MODELS_MIRROR_H_

#include "TmmBoundaryModel.h"
#include "Tmm.h"

namespace libMesh
{
  class Elem;
}
class TBDLLOCAL Mirror : public TmmBoundaryModel
{

  public:

    //! Destructor
    ~Mirror(void) {};

    //! Creator function
    static Mirror* create(const ModelOptions& options);


    virtual void Calculate_M_Matrix(void);

  protected:

    virtual void do_init(void);




  private:

    //! Constructor
    Mirror(const ModelOptions& options);
    double _member00;
    double _member01;
    double _member10;
    double _member11;


};



inline
Mirror::Mirror(const ModelOptions& options) :
  TmmBoundaryModel(options)
{
}



inline
Mirror*
Mirror::create(const ModelOptions& options)
{
  return new Mirror(options);
}


#endif /* SRC_CORE_MODULES_TMM_MODELS_MIRROR_H_ */
