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
 *   R_{aug} = C_nn(p^2-n_i^2) + C_pp(n^2-n_i^2)
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

    //! The parameter Cn
    double _Cn;

    //! The parameter Cp
    double _Cp;

};


//
// inline methods
// 


inline
AugerRecombination::AugerRecombination(void)
  : _Cn(1e-30),
    _Cp(1e-30)
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
  _Cn = mod->_Cn;
  _Cp = mod->_Cp;
}


#endif // _AUGERRECOMBINATION_H_
