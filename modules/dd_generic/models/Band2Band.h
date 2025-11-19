// $Id: Band2Band.h 3417 2012-09-12 16:16:07Z maufder $

#ifndef _BAND2BAND_H_
#define _BAND2BAND_H_

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"


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

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

    
  protected:

    //! Constructor
    Band2Band(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);


  private:

    //! B parameter in cm^-1/2 * V^-5/2 * s^-1
    double _B_param;

    //! Critical field in V/cm
    double _E0;

    //! Exponent for field dependency
    /*!
     * = 2 for direct transition, = 5/2 for indirect
     */
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
