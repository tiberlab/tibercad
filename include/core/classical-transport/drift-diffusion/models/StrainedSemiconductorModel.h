// $Id$

#ifndef _STRAINEDSEMICONDUCTOR_H_
#define _STRAINEDSEMICONDUCTOR_H_

#include "SemiconductorModel.h"

#include <map>

class Elem;
class Macrostrain;

//! A drift-diffusion model for a strained semiconductor
/*!
 * This implementation of a strained semiconductor is based on
 * the \c SemiconductorModel class. Based on a \c Macrostrain object
 * it will recalculate for every element the most important equilibrium
 * properties and consider strain induced polarization.
 */
class StrainedSemiconductorModel : public SemiconductorModel
{

  public:

    //! The default constructor
    StrainedSemiconductorModel(void);
    
    //! The destructor
    virtual ~StrainedSemiconductorModel(void) {};

    // Create a StrainedSemiconductorModel object
    static StrainedSemiconductorModel* create(void);

    //! Ignore strain related effects
    void ignore_strain(void);

    //! Include strain related effects
    void include_strain(void);

    //! Clean the internal cache of element data
    /*!
     * Band and equilibrium parameters are cached for each element so they
     * don't have to be recalculated during drift diffusion solving steps
     */
    void reset(void);


  protected:

    //! \copydoc DriftDiffusionProperties::prepare_element_data()
    virtual void prepare_element_data(void);

    //! \copydoc DriftDiffusionProperties::do_init()
    virtual void do_init();

    //! \copydoc DriftDiffusionProperties::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    //! \copydoc DriftDiffusionProperties::copy_from()
    virtual void copy_from(const PhysicalModelInterface* rhs);

  private:

    //! The data structure for the cached data
    struct ElementData
    {
      double Ec;
      double Ev;
      double mc;
      double mv;

      double Ef0;
      double ni;

      RealVectorValue polarization;
    };

    typedef std::map<const Elem*, ElementData> DataMap;
    
    StrainedSemiconductorModel(const StrainedSemiconductorModel& model);
    StrainedSemiconductorModel&
      operator=(const StrainedSemiconductorModel& model);
    
    Macrostrain* strain_model_;
    bool ignore_strain_;

    DataMap element_data_;

};



//
// inline methods
//

inline
PhysicalModelInterface*
StrainedSemiconductorModel::create_new(void) const
{
  return new StrainedSemiconductorModel();
}


inline
StrainedSemiconductorModel*
StrainedSemiconductorModel::create(void)
{
  return new StrainedSemiconductorModel();
}


inline
void
StrainedSemiconductorModel::ignore_strain(void)
{
  ignore_strain_ = true;
}


inline
void
StrainedSemiconductorModel::include_strain(void)
{
  ignore_strain_ = false;
}



#endif
