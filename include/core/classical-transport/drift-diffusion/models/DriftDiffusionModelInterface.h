// $Id$

#ifndef _DRIFTDIFFUSIONMODELINTERFACE_H_
#define _DRIFTDIFFUSIONMODELINTERFACE_H_

#include "TypeDefs.h"
#include "PhysicalModelInterface.h"

#include <cassert>
#include <map>

class DriftDiffusionProperties;

//! Base class for the different models used in Drift-Diffusion calculations
/*!
 * This is the base class for all implementations of mobility models,
 * recombination models etc. which will be used in conjunction with
 * a semiconductor model derived from DriftDiffusionProperties
 */
class DriftDiffusionModelInterface : public PhysicalModelInterface
{

  public:

    //! Destructor
    virtual ~DriftDiffusionModelInterface(void) {};

    //! Set the link to the DriftDiffusionProperties object
    /*!
     * \param dd_prop a pointer to the DriftDiffusionProperties object this
     * model belongs to
     */
    void set_driftdiffusionproperties(DriftDiffusionProperties* dd_prop);

    //! Get a reference to the DriftDiffusionProperties object
    /*!
     * \return a reference to the DriftDiffusionProperties object this
     * model belongs to
     */
    DriftDiffusionProperties& get_driftdiffusionproperties(void);


  protected:

    //! Empty constructor
    DriftDiffusionModelInterface(void);

    //! \copydoc PhysicalModelInterfaceInterface::copy_from()
    virtual void copy_from(const PhysicalModelInterface* rhs);

    //! The standard reference temperature in eV
    /*!
     * We give this value here because we get the lattice temperature
     * from DriftDiffusionProperties in eV and many models need something
     * like T/T0
     */
    static const double T0;


  private:

    //! Disable copy constructor
    DriftDiffusionModelInterface(const DriftDiffusionModelInterface&);

    //! Disable assignment operator
    DriftDiffusionModelInterface& operator=(const DriftDiffusionModelInterface&);
    
    //! The DriftDiffusionProperties object this model belongs to
    DriftDiffusionProperties* _dd_prop;

};


inline
DriftDiffusionModelInterface::DriftDiffusionModelInterface(void)
  : _dd_prop(0)
{
}


inline
void
DriftDiffusionModelInterface::set_driftdiffusionproperties(
    DriftDiffusionProperties* dd_prop)
{
  assert(dd_prop != 0);
  _dd_prop = dd_prop;
}


inline
DriftDiffusionProperties&
DriftDiffusionModelInterface::get_driftdiffusionproperties(void)
{
  assert(_dd_prop != 0);
  return *_dd_prop;
}

inline
void
DriftDiffusionModelInterface::copy_from(const PhysicalModelInterface* rhs)
{
  ignore_unused_variable(rhs);
}


#endif // _DRIFTDIFFUSIONMODELINTERFACE_H_
