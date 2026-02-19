/*  
 * This file is part of the tiberCAD module common.
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
 * \file DeltaDOS.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef TC_DELTADOS_H
#define TC_DELTADOS_H


#include "tibercad/physics/misc/DensityOfStates.h"


/*!
 * \brief Base class for density of states
 */
class TC_DLEXPORT DeltaDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~DeltaDOS(void) {};


    //! Creator function
    static DeltaDOS* create(const ModelOptions& options);




  protected:

    //! Constructor
    DeltaDOS(const ModelOptions& options);

    virtual void do_init(void);


    //! Get occupied states and derivative w.r.t. Ef
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;


  private:

  //! The effective DOS
  double _N0;
};

//
// inline methods
//

inline
DeltaDOS*
DeltaDOS::create(const ModelOptions& options)
{
  return new DeltaDOS(options);
}


#endif // TC_DELTADOS_H
