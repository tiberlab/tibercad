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
 * \file TTARecombination.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_TTARECOMBINATION_H
#define TC_TTARECOMBINATION_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of triplet-triplet annihilation
/*!
 * Triplet-triplet annihilation (TTA) is
 * modeled by \f[R=C n_T^{2}(1-exp{-2/kT\phi_T}\f]
 */
class TC_DLLOCAL TTARecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~TTARecombination(void) {};

    //! Create a ConstantMobility object
    static TTARecombination* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    TTARecombination(const ModelOptions& options);


    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;


    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Recombination rate parameter
    double C_;

};



//
// inline methods
// 

inline
TTARecombination::TTARecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0)
{
}


inline
TTARecombination*
TTARecombination::create(const ModelOptions& options)
{
  return new TTARecombination(options);
}






#endif // TC_TTARECOMBINATION_H
