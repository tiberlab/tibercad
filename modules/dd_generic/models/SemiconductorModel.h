// $Id: SemiconductorModel.h 3532 2013-02-20 12:34:19Z maufder $

#ifndef _SEMICONDUCTORMODEL_H_
#define _SEMICONDUCTORMODEL_H_

#include "tibercad/base/SimulationOptions.h"
#include "DDBulkModel.h"
#include "tibercad/physics/StrainInterface.h"

#include <vector>
#include <string>


// forward declarations
//class Point;
class Elem;


//! The default DD semiconductor model container
class TBDLLOCAL SemiconductorModel : public DDBulkModel
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

      libMesh::RealVectorValue polarization;
    };

    //! A data map type
    typedef std::map<const Elem*, ElementData> DataMap;
  

    /*! \copydoc DriftDiffusionProperties::do_init() */
    virtual void do_init();


    //virtual void prepare_element_data(void);


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
SemiconductorModel::DataMap&
SemiconductorModel::get_data_map(void)
{
  return _element_data;
}


#endif //_SEMICONDUCTORMODEL_H_
