// $Id: DriftDiffusionModelInterface.h 2882 2011-07-18 10:09:45Z maufder $

#ifndef _DRIFTDIFFUSIONMODELINTERFACE_H_
#define _DRIFTDIFFUSIONMODELINTERFACE_H_


#include "TypeDefs.h"
#include "PhysicalModel.h"

#include <cassert>

class DriftDiffusionProperties;

//! Base class for the different models used in Drift-Diffusion calculations
/*!
 * This is the base class for all implementations of mobility models,
 * recombination models etc. which will be used in conjunction with
 * a semiconductor model derived from DriftDiffusionProperties
 */
class TBDLEXPORT DriftDiffusionModelInterface : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~DriftDiffusionModelInterface(void) {};


    //! Get a reference to the DriftDiffusionProperties object
    /*!
     * \return a reference to the DriftDiffusionProperties object this
     * model belongs to
     */
    DriftDiffusionProperties& get_driftdiffusionproperties(void) const;


    //! Get a reference to the DriftDiffusionProperties object
    /*!
     * \return a reference to the DriftDiffusionProperties object this
     * model belongs to
     */
    DriftDiffusionProperties& get_bulk_driftdiffusionproperties(void) const;


  protected:

    //! Empty constructor
    DriftDiffusionModelInterface(const ModelOptions& options);

    //! The standard reference temperature in eV
    /*!
     * We give this value here because we get the lattice temperature
     * from DriftDiffusionProperties in eV and many models need something
     * like T/T0
     */
    static const double T0;


  private:

    
};


inline
DriftDiffusionModelInterface::DriftDiffusionModelInterface(const ModelOptions& options)
  : PhysicalModel(options)
{
}






#endif // _DRIFTDIFFUSIONMODELINTERFACE_H_
