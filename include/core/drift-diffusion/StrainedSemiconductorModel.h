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
    StrainedSemiconductorModel(Macrostrain* strain);
    
    //! The destructor
    virtual ~StrainedSemiconductorModel(void) {};

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

    virtual void prepare_element_data(void);

  private:

    //! The data structure for the cached data
    struct ElementData
    {
      double Ec;
      double Ev;
      double mc;
      double mv;

      double Ef0;
      double n0;
      double p0;

      RealVectorValue polarization;
    };

    typedef std::map<const Elem*, ElementData> DataMap;
    
    StrainedSemiconductorModel(const StrainedSemiconductorModel& model);
    
    Macrostrain* _strain;
    bool _ignore_strain;

    DataMap _element_data;

};

inline
void
StrainedSemiconductorModel::ignore_strain(void)
{
  _ignore_strain = true;
}

inline
void
StrainedSemiconductorModel::include_strain(void)
{
  _ignore_strain = false;
}


#endif
