// $Id$

#ifndef _EPQUENCHING_H_
#define _EPQUENCHING_H_

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of a triplet-polaron quenching
/*!
 * Triplet-polaron quenching (EPQ) is modeled as
 * \f[R=C n_i n_T (1-exp^{-\phi_T / kT}\f]
 * where \f$n_i\f$ is the quenching particle density (e.g. electron)
 */
class TBDLLOCAL EPQuenching : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~EPQuenching(void) {};

    //! Create a ConstantMobility object
    static EPQuenching* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    EPQuenching(const ModelOptions& options);


    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;


    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Recombination rate parameter
    double C_;

    //! The quenching particle
    int _quencher;

};



//
// inline methods
// 

inline
EPQuenching::EPQuenching(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0)
{
}


inline
EPQuenching*
EPQuenching::create(const ModelOptions& options)
{
  return new EPQuenching(options);
}






#endif // _EPQUENCHING_H__
