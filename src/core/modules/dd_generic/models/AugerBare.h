// $Id$

#ifndef _AUGERBARE_H_
#define _AUGERBARE_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


//! Implementation of bare Auger process
/*!
 * This class implements the bare Auger recombination processes
 * involving three particles, e.g.
 * \[2n + p \rightleftarrow n^*\]
 */
class TBDLLOCAL AugerBare : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~AugerBare(void) {};

    //! Create a ConstantMobility object
    static AugerBare* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    AugerBare(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

    
  private:

    double _rate_constant;

};


//
// inline methods
// 


inline
AugerBare::AugerBare(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _rate_constant(1e-30)
{
}


inline
AugerBare*
AugerBare::create(const ModelOptions& options)
{
  return new AugerBare(options);
}





#endif // _AUGERBARE_H_
