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
 * \file AvalancheGeneration.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_AVALANCHEGENERATION_H
#define TC_AVALANCHEGENERATION_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"


//! Implementation of Impact Ionization model
/*!
 * This class implements the impact ionization model processes according to
 *
 * \f{eqnarray*}
 * G_{II} & =& \sum_i\alpha_i|j_i| \\
 * \alpha_i & = & \gamma a_i e^{-\frac{\gamma b_i}{|E|}} \\
 * \f}
 *
 */
class TC_DLLOCAL AvalancheGeneration : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~AvalancheGeneration(void);

    //! Create a ConstantMobility object
    static AvalancheGeneration* create(const ModelOptions& options);


  protected:

    //! Constructor
    AvalancheGeneration(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;


    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;


  private:

    //! Parameter \c a
    std::vector<double> _a_param;

    //! Parameter \c b
    std::vector<double> _b_param;

    //! Phonon energy
    double _w0;

};


//
// inline methods
//





inline
AvalancheGeneration*
AvalancheGeneration::create(const ModelOptions& options)
{
  return new AvalancheGeneration(options);
}




#endif // TC_AVALANCHEGENERATION_H
