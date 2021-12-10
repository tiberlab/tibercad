/*
 * Mirror.h
 *
 *  Created on: 30 Sep 2021
 *      Author: pamiri
 */

#ifndef SRC_CORE_MODULES_TMM_MODELS_INCoherance_H_
#define SRC_CORE_MODULES_TMM_MODELS_INCoherance_H_

#include "TmmBulkModel.h"


class TBDLLOCAL incoherance : public TmmBulkModel
{

  public:

    //! Destructor
    ~incoherance(void) {};

    //! Creator function
    static incoherance* create(const ModelOptions& options);

    double get_incoherance_index() const;



  protected:

    virtual void do_init(void);




  private:

    //! Constructor
    incoherance(const ModelOptions& options);



};



inline
incoherance::incoherance(const ModelOptions& options) :
  TmmBulkModel(options)
{
}



inline
incoherance*
incoherance::create(const ModelOptions& options)
{
  return new incoherance(options);
}


#endif /* SRC_CORE_MODULES_TMM_MODELS_incoherance_H_ */
