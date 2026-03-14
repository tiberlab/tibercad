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
 * \file ExcitonGeneration.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_EXCITONGENERATION_H
#define TC_EXCITONGENERATION_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"

class ExcitonTransport;


//! Implementation of Exciton generation
/*!
 * This class implements Exciton generation processes that can be
 * modeled by \f[G_{x}=Cnp\f]
 */
class TC_DLLOCAL ExcitonGeneration : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonGeneration(void) {};

    //! \copydoc RecombinationModelInterface::get_net_recombination_rates()
    void get_net_recombination_rates(double& recomb_e, double& recomb_h);

    /*!
     * \copydoc
     * RecombinationModelInterface::get_net_recombination_rate_derivatives()
     */
    void get_net_recombination_rate_derivatives(
        std::vector<double>& recomb_e, std::vector<double>& recomb_h);

    //! Set the direct recombination parameters
    void set_parameters(double C);


  protected:

    //! Constructor
    ExcitonGeneration(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void);

    
  private:

    //! Recombination rate parameter
    double C_;

};


//
// inline methods
// 



inline
void
ExcitonGeneration::set_parameters(double C)
{
  C_ = C;
}



#endif // TC_EXCITONGENERATION_H
