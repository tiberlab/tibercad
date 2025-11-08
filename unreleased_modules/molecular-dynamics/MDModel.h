// $Id: MyPoissonModel.h 2356 2011-02-19 22:54:42Z maufder $

#ifndef _MDMODEL_H_
#define _MDMODEL_H_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"

class Elem;


//! This is the base class for the Poisson physical model
class MDModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~MDModel(void) {};

    //! Creator function
    static MDModel* create(const Material* mat,const ModelOptions& options);


    //! Get the relative permittivity
    //const RealTensor& get_permittivity(void) const;


    //! Get the total polarization
   // const RealVectorValue& get_polarization(void) const;


    //! Get the total charge
    //double get_charge_density(void) const;


    //! Set the relative permittivity
    //RealTensor& get_permittivity(void);


    //! Set the total polarization
   // RealVectorValue& get_polarization(void);


    //! Set the total charge
    //void set_charge_density(double charge_density);

    //! Reinit for element \c elem
   // void set_element(const Elem* elem);


    //! Set the current point
   // void set_point(const Point& point);


    //! Calculate everything
    void calculate(void);


  protected:

    //! Constructor
    MDModel(const ModelOptions& options);

    virtual void do_init(void);

    virtual void create_submodels(void);

    //! do the actual calculation
    virtual void do_calculate(void);

    //const Elem* get_element(void) const { return _elem; }
  
    //const Point& get_point(void) const { return _point; }
  

  private:


    //! The polarization models
    //std::vector<PolarizationModel* > _pm;

    //! The element we are currently working on
    //const Elem* _elem;

    //! The point we are currently using
    //Point _point;

    //! The relative permittivity
    //RealTensor _permittivity;

    //! The polarization
    //RealVectorValue _polarization;

    //! The charge density
    //double _charge;

    //! The permittivity model
    //PermittivityModel* _permittivity_model;

    //! The charge density model
    //MyChargeDensityModel* _charge_density;

    //! The constructor method

    static TiberModelObject*  _create(const ModelOptions& options, const void*);



    //! The destructor method
    static void _destroy(TiberModelObject* p);

};


inline
TiberModelObject*  MDModel::_create(const ModelOptions& options, const void*)
{
  return new MDModel(options);
}


inline
void
MDModel::_destroy(TiberModelObject* p)
{
  delete p;
}


//inline
//void
//MDModel::set_element(const Elem* elem)
//{
//  _elem = elem;
//}

//inline
//void
//MDModel::set_point(const Point& point)
//{
//  _point = point;
//}


inline
void
MDModel::calculate(void)
{
  do_calculate();
}

#endif // _MDMODEL_H_
