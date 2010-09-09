// $Id$

#ifndef _THERMALMODEL_H_
#define _THERMALMODEL_H_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "HeatTransportModel.h"
//#include "tiber_dll.h"


class Elem;
class HeatSourceModel;
class ThermalConductivityModel;

//class TBDLLOCAL ThermalModel : public PhysicalModel
class ThermalModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~ThermalModel(void) {};

    //! Creator function
    static ThermalModel* create(const ModelOptions& options);

  HeatTransportModel* get_heat_transport_model(void) const;

  //ThermalConductivityModel* get_thermal_conductivity_model(void) const;


  const double get_total_heat_source(void) const;

  const RealTensor& get_total_thermal_conductivity(void) const;

  const double get_sound_velocity(void) const;

  const double get_heat_capacity(void) const;

  const double get_relaxation_time(void) const;

  //! Reinit for element \c elem
  // void set_element(const Elem* elem);
  
  //! Set the current point
  //  void set_point(const Point& point);
  
  //! Calculate everything
  //void calculate(void);

 //! Calculate for a point on the given side
  virtual void calculate(const Elem* elem, const Point& point);

  protected:

    //! Constructor
    ThermalModel(const ModelOptions& options); 
    //! Create a new instance of this type
    virtual PhysicalModelInterface* create_new(void) const;


  void do_init_alloy(const PhysicalModelInterface *comp_A,
			    const PhysicalModelInterface *comp_B, double xa);
    //! Initialize
    virtual void do_init(void);

    //! Read database
    virtual void read_database(void);

//     //! do the actual calculation
//     virtual void do_calculate(void);

    virtual void create_submodels(void);

 
  void do_print_info(void);

  private:

  static TiberModelObject*  _create(const ModelOptions& options);

  static void  _destroy( TiberModelObject* p);

  //! The element we are currently working on
  // const Elem* _elem;
  
  //! The point we are currently using
  // Point _point;

  double  _heat_source;

  RealTensor _kappa;
  
  double _vg;

  double _cg;

  double _tg;

  // double _kg;

  HeatTransportModel*  _htm; 

  ThermalConductivityModel* _tcm; 

  //!Heat Source model map
  std::vector<HeatSourceModel* > _hsm;

};

inline
TiberModelObject*  ThermalModel::_create(const ModelOptions& options)
{
  return new ThermalModel(options);
}

inline
void  ThermalModel::_destroy( TiberModelObject* p)
{
  delete p;
}




inline
const double 
ThermalModel::get_total_heat_source(void) const
{

  return _heat_source;

}

inline
const RealTensor&
ThermalModel::get_total_thermal_conductivity(void) const
{
  return  _kappa;
}

inline
const double
ThermalModel::get_sound_velocity(void) const
{
  return  _vg;
}

inline
const double
ThermalModel::get_heat_capacity(void) const
{
  return  _cg;
}

inline
const double
ThermalModel::get_relaxation_time(void) const
{
  return  _tg;
}


inline
HeatTransportModel*
ThermalModel::get_heat_transport_model(void) const
{

  return _htm;

}


inline
ThermalModel*
ThermalModel::create(const ModelOptions& options)
{
  return dynamic_cast<ThermalModel*>(PhysicalModelInterface::create(_create,_destroy,options));
}

inline
PhysicalModelInterface*
ThermalModel::create_new(void) const
{
  return new ThermalModel(get_options());
}


#endif // _DEFAULTMODEL_H_
