// $Id$

#ifndef _DEFAULTMODEL_H_
#define _DEFAULTMODEL_H_

#include "MyPoissonModel.h"



class DefaultModel : public MyPoissonModel
{

  public:

    //! Destructor
    virtual ~DefaultModel(void) {};

    //! Creator function
    static DefaultModel* create(const ModelOptions& options);


  protected:

    //! Constructor
    DefaultModel(const ModelOptions& options) : MyPoissonModel(options) {}

    //! Create a new instance of this type
    virtual PhysicalModelInterface* create_new(void) const;


    //! Initialize
    virtual void do_init(void);

    //! Read database
    virtual void read_database(void);

    //! do the actual calculation
    virtual void do_calculate(void);


  private:


};


inline
DefaultModel*
DefaultModel::create(const ModelOptions& options)
{
  return new DefaultModel(options);
}


inline
PhysicalModelInterface*
DefaultModel::create_new(void) const
{
  return new DefaultModel(get_options());
}

#endif // _DEFAULTMODEL_H_
