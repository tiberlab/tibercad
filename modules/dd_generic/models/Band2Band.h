/*  
 * This file is part of the tiberCAD module dd_generic.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file Band2Band.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_BAND2BAND_H
#define TC_BAND2BAND_H

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







#endif // TC_DIRECTRECOMBINATION_H
