/*  
 * This file is part of the tiberCAD module thermal.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file ThermalModel.h
 * \brief tiberCAD thermal module header.
 *
 * \note This file is part of module thermal.
 */


#ifndef TC_THERMALMODEL_H
#define TC_THERMALMODEL_H

#include "tibercad/physics/PhysicalModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/base/libMeshDefs.h"

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




#endif // TC_DEFAULTMODEL_H
