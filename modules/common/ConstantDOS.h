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
 * \file ConstantDOS.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef TC_CONSTANTDOS_H
#define TC_CONSTANTDOS_H


#include "tibercad/physics/misc/DensityOfStates.h"

class TBDLEXPORT ConstantDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~ConstantDOS(void) {};

    //! Creator function
    static ConstantDOS* create(const ModelOptions& options);

  protected:

    //! Constructor
    ConstantDOS(const ModelOptions& options);

    virtual void read_database(void);

    //! Get occupied states and the derivative with respect to phi
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double E, double Epot,
                                      double kT, double kTlattice, const Elem* elem, const Point& p) const;

    virtual void do_init(void);
	
  private:
  	
    //
    double _Ewidth;

    //Total density parameter
    double _N0;
};

//
// inline methods
//

inline
ConstantDOS*
ConstantDOS::create(const ModelOptions& options)
{
  return new ConstantDOS(options);
}


#endif // TC_CONSTANTDOS_H
