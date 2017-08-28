// $Id: ExcitonDecay.h 3414 2012-09-10 20:40:28Z maufder $

#ifndef _EXCITONDECAY_H_
#define _EXCITONDECAY_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class SimulationInterface;

//! Implementation of direct recombination
/*!
 * This class implements a recombination processes that can be
 * modeled by \f[R_{decay}=Cn(1 - exp(E_{F,n}/k_BT))\f] where \c n is the density of the carriers
 */
class TBDLLOCAL ExcitonDecay : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonDecay(void) {};

    //! Create a ConstantMobility object
    static ExcitonDecay* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    ExcitonDecay(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;

  private:


    //! Exciton recombination time
    double  _tau;


};



//
// inline methods
// 

inline
ExcitonDecay::ExcitonDecay(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _tau(1e9)
{
}


inline
ExcitonDecay*
ExcitonDecay::create(const ModelOptions& options)
{
  return new ExcitonDecay(options);
}






#endif // _EXCITONDECAY_H__
