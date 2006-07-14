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

    //! Get the ID of the model with name \c model_name
    /*!
     * If the model was not already registered, it gets inserted in the
     * model list and assigned a new ID.
     *
     * \param model_name the name of the model
     * \return the ID of the model or 0 if it does not exist
     */
    static ID get_id(const std::string& model_name);

    //! Creates a new named model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     * 
     * \param name the model name
     * \return a pointer to the newly created object
     */
    static DriftDiffusionModelInterface* create(const std::string& name);

    //! Creates a new named model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     * 
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static DriftDiffusionModelInterface* create(const std::string& name,
        const ModelOptions& options);

    //! Set options for this model
    virtual void set_model_options(const ModelOptions& options) {};


  protected:

    //! Empty constructor
    DriftDiffusionModelInterface(void);


  private:

    typedef std::map<const std::string, ID>::iterator model_id_iterator;

    //! The unique ID of this model
    ID _id;

    //! A map with ID/model name pairs
    /*!
     * Models are counted starting from 1. 0 means undefined model.
     */
    static std::map<const std::string, ID> _model_ids;

    //! The DriftDiffusionProperties object this model belongs to
    const DriftDiffusionProperties* _dd_prop;

    //! Register a new model
    /*!
     * This method registers every new model that gets created and assigns
     * it a unique model ID.
     */
    static void register_model(DriftDiffusionModelInterface* model);
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
