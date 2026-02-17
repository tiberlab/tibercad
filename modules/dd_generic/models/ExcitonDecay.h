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
 * \file ExcitonDecay.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_EXCITONDECAY_H
#define TC_EXCITONDECAY_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of direct recombination
/*!
 * This class implements a recombination processes that can be
 * modeled by \f[R_{decay}=Cn(1 - exp(E_{F,n}/k_BT))\f] where \c n is the density of the carriers
 */
class TBDLLOCAL ExcitonDecay : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonDecay(void) {};

    //! Create a ConstantMobility object
    static ExcitonDecay* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    ExcitonDecay(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;

  private:


    //! Exciton recombination time
    double  _tau;


};



//
// inline methods
// 

inline
ExcitonDecay::ExcitonDecay(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _tau(1e9)
{
}


inline
ExcitonDecay*
ExcitonDecay::create(const ModelOptions& options)
{
  return new ExcitonDecay(options);
}






#endif // TC_EXCITONDECAY_H
