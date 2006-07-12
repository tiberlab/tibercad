// $Id$

#ifndef _RECOMBINATIONMODELINTERFACE_H_
#define _RECOMBINATIONMODELINTERFACE_H_

#include "DriftDiffusionModelInterface.h"

#include <vector>

class RecombinationModelInterface : public DriftDiffusionModelInterface
{

  public:

    //! Destructor
    virtual ~RecombinationModelInterface(void) {};

    //! Get the electron and hole net recombination rates
    /*!
     * The net recombination rate is defined as \f$R_{net}=R-G\f$
     */
    virtual void get_net_recombination_rates(double& recomb_e,
        double& recomb_h) = 0;

    //! Get the electron and hole net recombination rate derivatives
    /*!
     * The net recombination rate is defined as \f$R_{net}=R-G\f$
     */
    virtual void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h) = 0;

    //! Creates a new named recombination model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     * 
     * \param name the model name
     * \return a pointer to the newly created object
     */
    static RecombinationModelInterface* create(const std::string& name);

    //! Creates a new named recombination model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     * 
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static RecombinationModelInterface* create(const std::string& name,
        const ModelOptions& options);

    //! Set options for this model
    virtual void set_model_options(const ModelOptions& options) {};

  protected:

    //! Empty constructor
    RecombinationModelInterface(void);

  private:

};


#endif // _RECOMBINATIONMODELINTERFACE_H_
