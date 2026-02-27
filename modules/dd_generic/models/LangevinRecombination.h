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
 * \file LangevinRecombination.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_LANGEVINRECOMBINATION_H
#define TC_LANGEVINRECOMBINATION_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of Langevin recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{Langevin}=\gamma \frac{2}{\varepsilon_r \varepsilon_0}(\mu_e + \mu_h)np\f]
 */
class TC_DLLOCAL LangevinRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~LangevinRecombination(void) {};


  protected:

    //! Constructor
    LangevinRecombination(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual double calculate_rate_and_derivatives(std::vector<double>& dPotentials) override;

  private:

    typedef std::map<std::pair<SimulationInterface*, SimulationInterface*>,
        std::pair<unsigned int, double> > QRecMap;

    //! Relative permittivity from database
    double _er;

    //! gamma factor
    double _gamma;

};



//
// inline methods
// 

inline
LangevinRecombination::LangevinRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _er(0.0),
    _gamma(1.0)
{
}



#endif // TC_LANGEVINRECOMBINATION_H
