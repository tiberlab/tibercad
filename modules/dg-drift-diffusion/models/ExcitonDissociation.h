// $Id: ExcitonDissociation.h 2399 2011-02-23 15:43:57Z maufder $

#ifndef _EXCITONDISSOCIATION_H_
#define _EXCITONDISSOCIATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class SimulationInterface;


//! Implementation of Exciton dissociation
/*!
 * This class implements Exciton dissociation process
 */
class TBDLLOCAL ExcitonDissociation : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonDissociation(void) {};

    //! Create a ConstantMobility object
    static ExcitonDissociation* create(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    virtual void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    virtual void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    
  protected:

    //! Constructor
    ExcitonDissociation(const ModelOptions& options);
    
    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModel* create_new(void) const;


  private:

    //! Damping factor
    double d_;

    //! The \c ExcitonTransport to use
    SimulationInterface* _exciton_sim;

    //! The ID of the needed variable
    ID _Rdiss_id;

};



//
// inline methods
// 

inline
ExcitonDissociation::ExcitonDissociation(const ModelOptions& options)
  : RecombinationModelInterface(options),
    d_(1.0),
    _exciton_sim(NULL)
{
}

inline
ExcitonDissociation*
ExcitonDissociation::create(const ModelOptions& options)
{
  return new ExcitonDissociation(options);
}



inline
PhysicalModel*
ExcitonDissociation::create_new(void) const
{
  return new ExcitonDissociation(get_options());
}


#endif // _EXCITONDISSOCIATION_H_
