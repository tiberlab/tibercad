// $Id: BoltzmannModel.h 2457 2011-03-06 23:52:12Z gromano $

#ifndef _THERMALMODEL_H_
#define _THERMALMODEL_H_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "HeatTransportModel.h"
#include "tiber_dll.h"


class Elem;
class HeatSourceModel;
class ThermalConductivityModel;

class TBDLLOCAL BoltzmannModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~BoltzmannModel(void) {};

    //! Creator function
    static BoltzmannModel* create(const Material* mat, const ModelOptions& options);


  HeatTransportModel* get_heat_transport_model(void) const;

  //ThermalConductivityModel* get_thermal_conductivity_model(void) const;

  const double get_total_heat_source(void) const;

  const RealTensor& get_total_thermal_conductivity(void) const;

  const double get_sound_velocity(void) const;

  const double get_heat_capacity(void) const;

  const double get_relaxation_time(void) const;


 //! Calculate for a point on the given side
  void calculate(const Elem* elem, const Point& point);

  protected:

    //! Constructor
    BoltzmannModel(const ModelOptions& options); 
    //! Create a new instance of this type
    virtual PhysicalModelInterface* create_new(void) const;


  void do_init_alloy(const PhysicalModelInterface *comp_A,
			    const PhysicalModelInterface *comp_B, double xa);
    //! Initialize
    virtual void do_init(void);

    //! Read database
  virtual void read_database(void){};

    virtual void prepare_submodels(void);

 
  void do_print_info(void);

  private:

  static TiberModelObject* _create(const ModelOptions& options, const void*);

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
TiberModelObject*  BoltzmannModel::_create(const ModelOptions& options, const void*)
{
  return new BoltzmannModel(options);
}



inline
void  BoltzmannModel::_destroy( TiberModelObject* p)
{
  delete p;
}




inline
const double 
BoltzmannModel::get_total_heat_source(void) const
{

  return _heat_source;

}

inline
const RealTensor&
BoltzmannModel::get_total_thermal_conductivity(void) const
{
  return  _kappa;
}

inline
const double
BoltzmannModel::get_sound_velocity(void) const
{
  return  _vg;
}

inline
const double
BoltzmannModel::get_heat_capacity(void) const
{
  return  _cg;
}

inline
const double
BoltzmannModel::get_relaxation_time(void) const
{
  return  _tg;
}


inline
HeatTransportModel*
BoltzmannModel::get_heat_transport_model(void) const
{

  return _htm;

}




inline
PhysicalModelInterface*
BoltzmannModel::create_new(void) const
{
  return new BoltzmannModel(get_options());
}


#endif // _DEFAULTMODEL_H_
