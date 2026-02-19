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
 * \file EPQuenching.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_EPQUENCHING_H
#define TC_EPQUENCHING_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;

//! Implementation of a triplet-polaron quenching
/*!
 * Triplet-polaron quenching (EPQ) is modeled as
 * \f[R=C n_i n_T (1-exp^{-\phi_T / kT}\f]
 * where \f$n_i\f$ is the quenching particle density (e.g. electron)
 */
class TC_DLLOCAL EPQuenching : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~EPQuenching(void) {};

    //! Create a ConstantMobility object
    static EPQuenching* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    EPQuenching(const ModelOptions& options);


    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;


    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Recombination rate parameter
    double C_;

    //! The quenching particle
    int _quencher;

};



//
// inline methods
// 

inline
EPQuenching::EPQuenching(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0)
{
}


inline
EPQuenching*
EPQuenching::create(const ModelOptions& options)
{
  return new EPQuenching(options);
}






#endif // TC_EPQUENCHING_H
