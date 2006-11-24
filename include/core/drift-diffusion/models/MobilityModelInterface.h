// $Id$

#ifndef _MOBILITYMODELINTERFACE_H_
#define _MOBILITYMODELINTERFACE_H_

#include "DriftDiffusionModelInterface.h"

#include <vector>

class MobilityModelInterface : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~MobilityModelInterface(void);

    //! Get the mobility
    virtual double get_mobility(void) = 0;

    //! Get the derivatives of the mobility
    virtual void get_mobility_derivatives(std::vector<double>& dm) = 0;

    //! Creates a new named mobility model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     * 
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static MobilityModelInterface* create(const std::string& name,
        const ModelOptions& options = ModelOptions());


  protected:

    //! \copydoc DriftDiffusionProperties::DriftDiffusionProperties()
    MobilityModelInterface(void);

  private:

};



//
// inline methods
//


inline
MobilityModelInterface::MobilityModelInterface(void)
{
}


inline
MobilityModelInterface::~MobilityModelInterface(void)
{
}


inline
MobilityModelInterface*
MobilityModelInterface::create(const std::string& name,
    const ModelOptions& options)
{
  return dynamic_cast<MobilityModelInterface*>(
      PhysicalModelInterface::create("mob_" + name, options));
}



#endif // _MOBILITYMODELINTERFACE_H_
