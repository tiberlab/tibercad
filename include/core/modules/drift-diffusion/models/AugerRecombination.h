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
 *  C_{n,p} = (A + B\frac{T}{T_0} + C\left(\frac{T}{T_0})^2)
 *      (1 + H e^{-\{n,p\}/N_0})
 * \f]
 */
class TBDLEXPORT AugerRecombination : public RecombinationModelInterface
{

  public:

    //! Constructor
    TBDLLOCAL AugerRecombination(void);

    //! Destructor
    virtual ~AugerRecombination(void) {};

    //! Create a ConstantMobility object
    static AugerRecombination* create(void);

    /*! \copydoc RecombinationModelInterface::get_net_recombination_rates() */
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*! \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    
  protected:

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    //! \copydoc RecombinationModelInterface::copy_from()
    virtual void copy_from(const PhysicalModelInterface* rhs);

    /*! \copydoc RecombinationModelInterface::calculate_VCA() */
    virtual void calculate_VCA(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

    
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
AugerRecombination::AugerRecombination(void)
  : _An(6.7000e-32),
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
AugerRecombination::create(void)
{
  return new AugerRecombination();
}


inline
PhysicalModelInterface*
AugerRecombination::create_new(void) const
{
  return new AugerRecombination();
}


inline
void
AugerRecombination::copy_from(const PhysicalModelInterface* rhs)
{
  RecombinationModelInterface::copy_from(rhs);
  
  const AugerRecombination* mod = dynamic_cast<const AugerRecombination*>(rhs);
  _An = mod->_An;
  _Ap = mod->_Ap;
  _Bn = mod->_Bn;
  _Bp = mod->_Bp;
  _Cn = mod->_Cn;
  _Cp = mod->_Cp;
  _Hn = mod->_Hn;
  _Hp = mod->_Hp;
  _N0n = mod->_N0n;
  _N0p = mod->_N0p;
}


#endif // _AUGERRECOMBINATION_H_
