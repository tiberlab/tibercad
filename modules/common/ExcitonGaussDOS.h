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
 * \file ExcitonGaussDOS.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */

#ifndef TC_EXCITONGAUSSDOS_H
#define TC_EXCITONGAUSSDOS_H


#include "tibercad/physics/misc/DensityOfStates.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT ExcitonGaussDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~ExcitonGaussDOS(void) {};


    //! Creator function
    static ExcitonGaussDOS* create(const ModelOptions& options);



  protected:

    //! Constructor
    ExcitonGaussDOS(const ModelOptions& options);

    //! Read parameters from database
    virtual void read_database(void);

    virtual void do_init(void);

    virtual void do_reinit(void);

    virtual void do_print_info(void);


    //! Get occupied states and derivative w.r.t. Ef
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;



  private:

    //! The exciton spin (for convenience)
    double _J;

    //! The exciton binding energy (eV)
    double _R;

    //! The gaussian DOS variance
    double _sigma;

    //! The exciton energy without considering binding energy (= gap)
    double _energy;

    //! Trapezoidal integration
    double _trapez(std::vector<double>& x, std::vector<double>& y) const;

    //! Bose - Einstein occupation function
    double _f(double x, double Ef, double kT) const;

    double _order;


};

//
// inline methods
//

inline
ExcitonGaussDOS*
ExcitonGaussDOS::create(const ModelOptions& options)
{
  return new ExcitonGaussDOS(options);
}


#endif // TC_BULKDOS_H
