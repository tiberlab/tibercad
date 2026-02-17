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
 * \file ExcitonGeneration.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_EXCITONGENERATION_H
#define TC_EXCITONGENERATION_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of exciton generation/dissociation from free carrier gas
/*!
 * This class implements generation and dissociation of excitons
 * modeled by \f[R=\gamma np(1 - exp(\frac{E_{Fx} - E_{Fn} + E_{Fp}}{k_BT}))\f]
 *
 * In the input file, the carriers have to be provided in the order
 * \c electron, \c hole, \c exciton
 */
class TBDLLOCAL ExcitonGeneration : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonGeneration(void) {};

    //! Create a ConstantMobility object
    static ExcitonGeneration* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    ExcitonGeneration(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Recombination rate parameters
    double  _gamma;

    bool _stat_fac;

};



//
// inline methods
// 

inline
ExcitonGeneration::ExcitonGeneration(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _stat_fac(true)
{
}


inline
ExcitonGeneration*
ExcitonGeneration::create(const ModelOptions& options)
{
  return new ExcitonGeneration(options);
}






#endif // TC_EXCITONGENERATION_H
