/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file ExcitonDissociation.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_EXCITONDISSOCIATION_H
#define TC_EXCITONDISSOCIATION_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class SimulationInterface;


//! Implementation of Exciton dissociation
/*!
 * This class implements Exciton dissociation process
 */
class TC_DLLOCAL ExcitonDissociation : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonDissociation(void) {};

    //! Create a ConstantMobility object
    static ExcitonDissociation* create(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    virtual void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    virtual void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    
  protected:

    //! Constructor
    ExcitonDissociation(const ModelOptions& options);
    
    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc RecombinationModelInterface::create_new()
    virtual PhysicalModel* create_new(void) const;


  private:

    //! Damping factor
    double d_;

    //! The \c ExcitonTransport to use
    SimulationInterface* _exciton_sim;

    //! The ID of the needed variable
    ID _Rdiss_id;

};



//
// inline methods
// 


inline
ExcitonDissociation*
ExcitonDissociation::create(const ModelOptions& options)
{
  return new ExcitonDissociation(options);
}



inline
PhysicalModel*
ExcitonDissociation::create_new(void) const
{
  return new ExcitonDissociation(get_options());
}


#endif // TC_EXCITONDISSOCIATION_H
