// $Id$

#ifndef _TTARECOMBINATION_H_
#define _TTARECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"

class SimulationInterface;

//! Implementation of a generic recombination model with constant parameter
/*!
 * This class implements recombination processes that can be
 * modeled by \f[R=C\prod_i n_i^{\alpha_i}(1-1/kT\sum_i\pm\alpha_i\phi_i\f]
 */
class TBDLLOCAL TTARecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~TTARecombination(void) {};

    //! Create a ConstantMobility object
    static TTARecombination* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    TTARecombination(const ModelOptions& options);


    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;


    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Recombination rate parameter
    double C_;

};



//
// inline methods
// 

inline
TTARecombination::TTARecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0)
{
}


inline
TTARecombination*
TTARecombination::create(const ModelOptions& options)
{
  return new TTARecombination(options);
}






#endif // _TTARECOMBINATION_H__
