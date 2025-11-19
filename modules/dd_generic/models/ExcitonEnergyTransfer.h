

#ifndef _EXCITONENERGYTRANSFER_H_
#define _EXCITONENERGYTRANSFER_H_

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"


//! Implementation of direct recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{direct}=C(np-n_i^2)\f]
 */
class TBDLLOCAL ExcitonEnergyTransfer : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonEnergyTransfer(void) {};

    //! Create a ConstantMobility object
    static ExcitonEnergyTransfer* create(const ModelOptions& options);


  protected:

    //! Constructor
    ExcitonEnergyTransfer(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::do_reinit()
    virtual void do_reinit(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Donor material permittivity
    double _er;

    //! Donor exciton effective mass
    double _m;

    //! Forster radius
    double _Rf;

    //! Dexter radius
    double _Rd;

    //! Average donor-acceptor distance
    double _R_da;

    //! Total donor exciton lifetime
    double _tau;

    //! Donor exciton radiative lifetime
    double _tau_rad;

    //! Forster rate
    double _Kf;

    //! Dexter rate
    double _Kd;

    // donor and acceptor ids
    ID _id_d;
    ID _id_a;


};



//
// inline methods
// 

inline
ExcitonEnergyTransfer::ExcitonEnergyTransfer(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _Rf(3e-7),
    _Rd(1e-7),
    _R_da(1.5e-7),
    _er(1.0),
    _m(2.0),
    _tau(1e-9),
    _tau_rad(1e-9)
{
}


inline
ExcitonEnergyTransfer*
ExcitonEnergyTransfer::create(const ModelOptions& options)
{
  return new ExcitonEnergyTransfer(options);
}






#endif // _EXCITONENERGYTRANSFER_H__
