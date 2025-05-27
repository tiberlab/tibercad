// $Id$

#ifndef _POISSONMODEL_H_
#define _POISSONMODEL_H_

#include "PhysicalModel.h"


//! This is the base class for the WI physical model
class WIModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~WIModel(void) {};

    //! Creator function
    static WIModel* create(const Material* mat, const ModelOptions& options);

    //! Calculate everything
    void calculate(const Elem* elem, const Point& point);

    //! Get the water solubility
    double get_solubility(void) const;

    //! Get the water diffusivity
    double get_diffusivity(void) const;


  protected:

    //! Constructor
    WIModel(const ModelOptions& options);

    virtual void do_init(void) override;

    virtual void prepare_submodels(void) override;

  

  private:
 
    //! The water solubility in g/cm^3/P
    double _solubility = 0.45e-6;;
 
    //! The water diffusivity in cm^2/s
    double _diffusivity = 3.43e-7;

    //! The constructor method
    static TiberModelObject* _create(const ModelOptions& options, const void*);

    //! The destructor method
    static void _destroy(TiberModelObject* p);

};




inline
WIModel::WIModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


inline
double
WIModel::get_solubility(void) const
{
  return _solubility;
}

inline
double
WIModel::get_diffusivity(void) const
{
  return _diffusivity;
}


#endif // _POISSONMODEL_H_
