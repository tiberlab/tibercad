// $Id: ThermalModel.h 2457 2011-03-06 23:52:12Z gromano $

#ifndef _THERMALMODEL_H_
#define _THERMALMODEL_H_

#include "PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tiber_dll.h"
#include "libMeshDefs.h"

//class Elem;
class HeatSourceModel;
class ThermalConductivityModel;

class TBDLLOCAL ThermalModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~ThermalModel(void) {};

    //! Creator function
    static ThermalModel* create(const Material* mat, const ModelOptions& options);

  //ThermalConductivityModel* get_thermal_conductivity_model(void) const;

  double get_total_heat_source(void) const;

  const libMesh::RealTensor& get_total_thermal_conductivity(void) const;

 //! Calculate for a point on the given side
  void calculate(const Elem* elem, const Point& point, double temperature);

  protected:

    //! Constructor
    ThermalModel(const ModelOptions& options); 


    //! Initialize
    virtual void do_init(void);

    //! Read database
  virtual void read_database(void){};

    virtual void prepare_submodels(void);

 
  void do_print_info(void);

  private:

  static TiberModelObject*  _create(const ModelOptions& options, const void*);

  static void  _destroy( TiberModelObject* p);



  ThermalConductivityModel* _tcm;

  libMesh::RealTensor _kappa;

  double  _heat_source;

  //!Heat Source model map
  std::vector<HeatSourceModel* > _hsm;

};

inline
TiberModelObject*  ThermalModel::_create(const ModelOptions& options, const void*)
{
  return new ThermalModel(options);
}

inline
void  ThermalModel::_destroy( TiberModelObject* p)
{
  delete p;
}




inline
double
ThermalModel::get_total_heat_source(void) const
{
  return _heat_source;

}

inline
const libMesh::RealTensor&
ThermalModel::get_total_thermal_conductivity(void) const
{
  return  _kappa;
}




#endif // _DEFAULTMODEL_H_
