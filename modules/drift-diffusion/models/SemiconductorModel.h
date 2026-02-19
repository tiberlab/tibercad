/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file SemiconductorModel.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_SEMICONDUCTORMODEL_H
#define TC_SEMICONDUCTORMODEL_H

#include "tibercad/base/SimulationOptions.h"
#include "DDBulkModel.h"
#include "tibercad/physics/StrainInterface.h"

#include <vector>
#include <string>


// forward declarations
//class Point;
//class Elem;


//! The default DD semiconductor model container
class TC_DLLOCAL SemiconductorModel : public DDBulkModel
{

  public:
    
    //! The constructor
    SemiconductorModel(const ModelOptions& options);
    
    //! The destructor
    virtual ~SemiconductorModel(void);

    //! This method creates a SimpleSemiconductorModel object
    static SemiconductorModel* create(const ModelOptions& options);


    //! Clean the internal cache of element data
    /*!
     * Band and equilibrium parameters are cached for each element so they
     * don't have to be recalculated during drift diffusion solving steps
     */
    void reset(void);

    

  protected:


    //! The data structure for element-wise cached data
    struct ElementData
    {
      double Ec;
      double Ev;
      double mc;
      double mv;

      double Ef0;
      //double ni;

libMesh::RealVectorValue polarization;
    };

    //! A data map type
    typedef std::map<const libMesh::Elem*, ElementData> DataMap;
  

    /*! \copydoc DriftDiffusionProperties::do_init() */
    virtual void do_init();


    /*! \copydoc DriftDiffusionProperties::prepare_element_data() */
    virtual void prepare_element_data(void);


    //! Set the object to unprepared state
    void set_to_unprepared(void);


    //! Get the data map with the element wise cached data
    DataMap& get_data_map(void);


  private:

    typedef DDBulkModel Parent;
    
    SemiconductorModel(const SemiconductorModel& model);
    SemiconductorModel& operator=(const SemiconductorModel& model);


    //! A flag to tell the state of this object
    /*!
     * \c true means that all data is prepared and ready for use
     */
    bool _is_prepared;


    //! The map with the element wise data
    DataMap _element_data;


    //! Should we always recompute band parameters?
    /*!
     * Use this for selfconsistent simulations
     */
    bool _recompute_band_parameters;


};


//
// inline member functions
//


inline
SemiconductorModel*
SemiconductorModel::create(const ModelOptions& options)
{
  return new SemiconductorModel(options);
}


inline
void
SemiconductorModel::set_to_unprepared(void)
{
  _is_prepared = false;
}




inline
SemiconductorModel::DataMap&
SemiconductorModel::get_data_map(void)
{
  return _element_data;
}


#endif //_SEMICONDUCTORMODEL_H_
