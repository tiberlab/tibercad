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
 * \file AugerBare.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_AUGERBARE_H
#define TC_AUGERBARE_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"


//! Implementation of bare Auger process
/*!
 * This class implements the bare Auger recombination processes
 * involving three particles, e.g.
 * \[2n + p \rightleftarrow n^*\]
 */
class TC_DLLOCAL AugerBare : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~AugerBare(void) {};

    //! Create a ConstantMobility object
    static AugerBare* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    AugerBare(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

    
  private:

    double _rate_constant;

};


//
// inline methods
// 


inline
AugerBare::AugerBare(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _rate_constant(1e-30)
{
}


inline
AugerBare*
AugerBare::create(const ModelOptions& options)
{
  return new AugerBare(options);
}





#endif // TC_AUGERBARE_H
