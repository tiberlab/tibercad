// $Id$

#ifndef _BAND2BAND_H_
#define _BAND2BAND_H_

#include "RecombinationModelInterface.h"
#include "TypeDefs.h"


//! Implementation of band-to-band tunneling
/*!
 * This class implements band-to-band tunneling in a local
 * approximation
 */
class TBDLLOCAL Band2Band : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~Band2Band(void) {};

    //! Create a ConstantMobility object
    static Band2Band* create(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);


    
  protected:

    //! Constructor
    Band2Band(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);


  private:

    //! B parameter
    double _B_param;

    //! Critical field
    double _E0;

    //! Exponent for field dependency
    double _sigma;


};



//
// inline methods
// 

inline
Band2Band::Band2Band(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _B_param(0.0),
    _E0(1e7),
    _sigma(2.5)
{
}


inline
Band2Band*
Band2Band::create(const ModelOptions& options)
{
  return new Band2Band(options);
}







#endif // _DIRECTRECOMBINATION_H__
