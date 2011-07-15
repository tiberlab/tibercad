// $Id$

#ifndef _SEMICONDUCTORMODEL_H_
#define _SEMICONDUCTORMODEL_H_

#include "SimulationOptions.h"
#include "DriftDiffusionProperties.h"
#include "StrainInterface.h"

#include <vector>
#include <string>


// forward declarations
//class Point;
class Elem;
class DDsemiconductor;


//! A semiconductor model using k.p theory
/*!
 * This model calculates band properties using the \c DDsemiconductor
 * interface, which is based on k.p theory.
 * 
 */
class TBDLLOCAL SemiconductorModel : public DriftDiffusionProperties
{

  public:
    
    //! The constructor
    SemiconductorModel(const ModelOptions& options);
    
    //! The destructor
    virtual ~SemiconductorModel(void);

    //! This method creates a SimpleSemiconductorModel object
    static SemiconductorModel* create(const ModelOptions& options);

    /*! \copydoc DriftDiffusionProperties::calculate_equilibrium_properties() */
    virtual void calculate_equilibrium_properties(void);


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

      RealVectorValue polarization;
    };

    //! A data map type
    typedef std::map<const Elem*, ElementData> DataMap;
  

    /*! \copydoc DriftDiffusionProperties::do_init() */
    virtual void do_init();

    //! Read the from database
    virtual void create_submodels(void);


    //! Get the physical semiconductor model
    /*!
     * Derived classes will need to access the physical model, e.g. to
     * set the strain.
     */
    DDsemiconductor* get_physical_model(void)
      { return _bulk_model; };

    //! Extract the band properties from _bulk_model
    /*!
     * This method looks for the band extrema and puts the effective
     * mass, band edges etc. into the BandProperties structure
     */
    void extract_band_properties(void);

    /*! \copydoc DriftDiffusionProperties::prepare_element_data() */
    virtual void prepare_element_data(void);


    //! Set the object to unprepared state
    void set_to_unprepared(void);


    /*! \copydoc PhysicalModelInterface::do_print_info() */
    //virtual void do_print_info(void);


    //! Get the data map with the element wise cached data
    DataMap& get_data_map(void);


  private:

    typedef DriftDiffusionProperties Parent;
    
    SemiconductorModel(const SemiconductorModel& model);
    SemiconductorModel& operator=(const SemiconductorModel& model);


    //! A flag to tell the state of this object
    /*!
     * \c true means that all data is prepared and ready for use
     */
    bool _is_prepared;


    //! The physical model for this semiconductor
    /*!
     * The physical model is based on an effective mass approximation
     */
    DDsemiconductor* _bulk_model;


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
