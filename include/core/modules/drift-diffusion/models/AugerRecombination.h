// $Id$

#ifndef _AUGERRECOMBINATION_H_
#define _AUGERRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


//! Implementation of Auger recombination
/*!
 * This class implements Auger recombination processes that can be
 * modeled by 
 * \f[
 *   R_{aug} = (C_nn + C_pp)(np-n_i^2)
 * \f]
 * with
 * \f[
 *  C_{\{n,p\}} = \left(A + B\frac{T}{T_0} + C\left(\frac{T}{T_0}\right)^2\right)
 *      \left(1 + H e^{-\{n,p\}/N_0}\right)
 * \f]
 */
class TBDLLOCAL AugerRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~AugerRecombination(void) {};

    //! Create a ConstantMobility object
    static AugerRecombination* create(const ModelOptions& options);

    /*! \copydoc RecombinationModelInterface::get_net_recombination_rates() */
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*! \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    
  protected:

    //! Constructor
    AugerRecombination(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    
  private:

    //! The parameter A for the electrons
    double _An;

    //! The parameter A for the holes
    double _Ap;

    //! The parameter B for the electrons
    double _Bn;

    //! The parameter B for the holes
    double _Bp;

    //! The parameter C for the electrons
    double _Cn;

    //! The parameter C for the holes
    double _Cp;

    //! The parameter H for the electrons
    double _Hn;

    //! The parameter H for the holes
    double _Hp;

    //! The parameter N0 for the electrons
    double _N0n;

    //! The parameter N0 for the holes
    double _N0p;


    //! Whether to use fixed Cn or not
    bool _fixed_Cn;

    //! Whether to use fixed Cp or not
    bool _fixed_Cp;


    //! Get Cn
    double get_Cn(void);

    //! Get Cp
    double get_Cp(void);

};


//
// inline methods
// 


inline
AugerRecombination::AugerRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _An(6.7000e-32),
    _Ap(7.2000e-32),
    _Bn(2.4500e-31),
    _Bp(4.5000e-33),
    _Cn(-2.2000e-32),
    _Cp(2.6300e-32),
    _Hn(3.46667),
    _Hp(8.25688),
    _N0n(1.0000e+18),
    _N0p(1.0000e+18),
    _fixed_Cn(false),
    _fixed_Cp(false)
{
}


inline
AugerRecombination*
AugerRecombination::create(const ModelOptions& options)
{
  return new AugerRecombination(options);
}





#endif // _AUGERRECOMBINATION_H_
