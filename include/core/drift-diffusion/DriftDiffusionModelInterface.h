// $Id$

#ifndef _DRIFTDIFFUSIONMODELINTERFACE_H_
#define _DRIFTDIFFUSIONMODELINTERFACE_H_

#include "TypeDefs.h"

#include <cassert>
#include <map>

class DriftDiffusionProperties;

//! Base class for the different models used in Drift-Diffusion calculations
/*!
 * This is the base class for all implementations of mobility models,
 * recombination models etc. which will be used in conjunction with
 * a semiconductor model derived from DriftDiffusionProperties
 */
class DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~DriftDiffusionModelInterface(void) {};

    //! Set the link to the DriftDiffusionProperties object
    /*!
     * \param dd_prop a pointer to the DriftDiffusionProperties object this
     * model belongs to
     */
    void set_driftdiffusionproperties(const DriftDiffusionProperties* dd_prop);

    //! Get a reference to the DriftDiffusionProperties object
    /*!
     * \return a reference to the DriftDiffusionProperties object this
     * model belongs to
     */
    const DriftDiffusionProperties& get_driftdiffusionproperties(void) const;

    //! Get the unique ID of this model
    ID get_id(void) const;


  protected:

    //! Empty constructor
    DriftDiffusionModelInterface(void);

    //! Register a new model
    /*!
     * This method registers every new model that gets created and assigns
     * it a unique model ID.
     * This method has to be called from the constructor of any model that
     * can be instantiated
     */
    void register_model(void);

  private:

    typedef std::map<const std::string, ID>::iterator model_id_iterator;

    //! The unique ID of this model
    ID _id;

    //! A map with ID/model name pairs
    static std::map<const std::string, ID> _model_ids;

    //! The DriftDiffusionProperties object this model belongs to
    const DriftDiffusionProperties* _dd_prop;
};


inline
DriftDiffusionModelInterface::DriftDiffusionModelInterface(void)
  : _dd_prop(0)
{
}

inline
ID
DriftDiffusionModelInterface::get_id(void) const
{
  return _id;
}

inline
void
DriftDiffusionModelInterface::set_driftdiffusionproperties(
    const DriftDiffusionProperties* dd_prop)
{
  assert(dd_prop != 0);
  _dd_prop = dd_prop;
}

inline
const DriftDiffusionProperties&
DriftDiffusionModelInterface::get_driftdiffusionproperties(void) const
{
  assert(_dd_prop != 0);
  return *_dd_prop;
}



#endif // _DRIFTDIFFUSIONMODELINTERFACE_H_
