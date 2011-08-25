// $Id: ThermalModel.h 2457 2011-03-06 23:52:12Z gromano $

#ifndef _THERMALMODEL_H_
#define _THERMALMODEL_H_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"


class Elem;
class HeatSourceModel;
class ThermalConductivityModel;

class TBDLLOCAL ThermalModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~ThermalModel(void) {};

    //! Creator function
    static ThermalModel* create(const ModelOptions& options);


  //ThermalConductivityModel* get_thermal_conductivity_model(void) const;


  const double get_total_heat_source(void) const;

  const RealTensor& get_total_thermal_conductivity(void) const;

 //! Calculate for a point on the given side
  void calculate(const Elem* elem, const Point& point);

  protected:

    //! Constructor
    ThermalModel(const ModelOptions& options); 
    //! Create a new instance of this type
    virtual PhysicalModelInterface* create_new(void) const;


  void do_init_alloy(const PhysicalModelInterface *comp_A,
		     const PhysicalModelInterface *comp_B, double xa){};
    //! Initialize
    virtual void do_init(void);

    //! Read database
  virtual void read_database(void){};

    virtual void prepare_submodels(void);

 
  void do_print_info(void);

  private:

  static TiberModelObject*  _create(const ModelOptions& options);

  static void  _destroy( TiberModelObject* p);



  double  _heat_source;

  RealTensor _kappa;
  
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
