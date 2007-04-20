// $Id$

#ifndef _SEMICONDUCTORMODEL_H_
#define _SEMICONDUCTORMODEL_H_

#include "SimulationOptions.h"
#include "DriftDiffusionProperties.h"

#include <vector>
#include <string>


// forward declarations
//class Point;
class Elem;
class DDsemiconductor;

//! A generic semiconductor model
/*!
 * This model calculates band properties using the \c DDsemiconductor
 * interface. It is not intended for use in strained structures.
 * 
 */
class SemiconductorModel : public DriftDiffusionProperties
{
  public:
    
    //! The constructor
    /*!
     * This constructor can be called be derived classes to specify
     * a different model name.
     */
    SemiconductorModel(void);

    //! The destructor
    virtual ~SemiconductorModel(void);

    //! This method creates a SimpleSemiconductorModel object
    static SemiconductorModel* create(void);

    /*! \copydoc DriftDiffusionProperties::calculate_equilibrium_properties() */
    virtual void calculate_equilibrium_properties(void);

    /*! \copydoc DriftDiffusionProperties::calculate_VCA() */
    virtual void calculate_VCA(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

    //! \deprecated { Was for debugging }
    void print_info(void) const;

    
  protected:

    /*! \copydoc DriftDiffusionProperties::do_init() */
    virtual void do_init();

    //! Read the from database
    virtual void read_database(void);

    /*! \copydoc DriftDiffusionProperties::create_new() */
    virtual PhysicalModelInterface* create_new(void) const;

    /*! \copydoc DriftDiffusionProperties::copy_from() */
    virtual void copy_from(const PhysicalModelInterface* rhs);

    //! Get the physical semiconductor model
    /*!
     * Derived classes will need to access the physical model, e.g. to
     * set the strain.
     */
    DDsemiconductor* get_physical_model(void)
      { return bulk_model_; };

    //! Extract the band properties from bulk_model_
    /*!
     * This method looks for the band extrema and puts the effective
     * mass, band edges etc. into the BandProperties structure
     */
    void extract_band_properties(void);

    /*! \copydoc DriftDiffusionProperties::prepare_element_data() */
    virtual void prepare_element_data(void);


    //! Set the object to unprepared state
    void set_to_unprepared(void);

  private:

    typedef DriftDiffusionProperties Parent;
    
    SemiconductorModel(const SemiconductorModel& model);
    SemiconductorModel& operator=(const SemiconductorModel& model);

    //! A flag to tell the state of this object
    /*!
     * \c true means that all data is prepared and ready for use
     */
    bool is_prepared_;

    //! The physical model for this semiconductor
    /*!
     * The physical model is based on an effective mass approximation
     */
    DDsemiconductor* bulk_model_;


};


//
// inline member functions
//


inline
SemiconductorModel*
SemiconductorModel::create(void)
{
  return new SemiconductorModel();
}


inline
void
SemiconductorModel::set_to_unprepared(void)
{
  is_prepared_ = false;
}

inline
PhysicalModelInterface*
SemiconductorModel::create_new(void) const
{
  return new SemiconductorModel();
}


inline
void
SemiconductorModel::copy_from(const PhysicalModelInterface* rhs)
{
  Parent::copy_from(rhs);

  const SemiconductorModel* mod = dynamic_cast<const SemiconductorModel*>(rhs);
  is_prepared_ = mod->is_prepared_;

}


#endif //_SEMICONDUCTORMODEL_H_
