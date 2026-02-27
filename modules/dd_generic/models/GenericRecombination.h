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
 * \file GenericRecombination.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_GENERICRECOMBINATION_H
#define TC_GENERICRECOMBINATION_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of a generic recombination model with constant parameter
/*!
 * This class implements recombination processes that can be
 * modeled by \f[R=C\prod_i n_i^{\alpha_i}(1-\exp{1/kT\sum_i\pm\alpha_i\phi_i}\f]
 */
class TC_DLLOCAL GenericRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~GenericRecombination(void) {};

    
  protected:

    //! Constructor
    GenericRecombination(const ModelOptions& options);


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
GenericRecombination::GenericRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0)
{
}






#endif // TC_GENERICRECOMBINATION_H
