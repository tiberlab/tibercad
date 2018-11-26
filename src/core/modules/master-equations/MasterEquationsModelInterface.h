// $Id: DriftDiffusionModelInterface.h 2882 2011-07-18 10:09:45Z maufder $

#ifndef _MASTEREQUATIONSMODELINTERFACE_H_
#define _MASTEREQUATIONSMODELINTERFACE_H_


#include "TypeDefs.h"
#include "PhysicalModelInterface.h"

#include <cassert>

class MasterEquationsProperties;

//! Base class for the different models used in Drift-Diffusion calculations
/*!
 * This is the base class for all implementations of mobility models,
 * recombination models etc. which will be used in conjunction with
 * a semiconductor model derived from DriftDiffusionProperties
 */
class TBDLEXPORT MasterEquationsModelInterface : public PhysicalModelInterface
{

  public:

    //! Destructor
    virtual ~MasterEquationsModelInterface(void) {};


    //! Get a reference to the MasterEquationsProperties object
    /*!
     * \return a reference to the MasterEquationsProperties object this
     * model belongs to
     */
    MasterEquationsProperties& get_masterequationsproperties(void) const;


  protected:

    //! Empty constructor
    MasterEquationsModelInterface(const ModelOptions& options);

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
MasterEquationsModelInterface::MasterEquationsModelInterface(const ModelOptions& options)
  : PhysicalModelInterface(options)
{
}






#endif // _MASTEREQUATIONSMODELINTERFACE_H_
