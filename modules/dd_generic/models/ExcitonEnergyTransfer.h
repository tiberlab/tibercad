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
 * \file ExcitonEnergyTransfer.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */



#ifndef TC_EXCITONENERGYTRANSFER_H
#define TC_EXCITONENERGYTRANSFER_H

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"


//! Implementation of direct recombination
/*!
 * This class implements direct recombination processes that can be
 * modeled by \f[R_{direct}=C(np-n_i^2)\f]
 */
class TC_DLLOCAL ExcitonEnergyTransfer : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~ExcitonEnergyTransfer(void) {};


  protected:

    //! Constructor
    ExcitonEnergyTransfer(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::do_reinit()
    virtual void do_reinit(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials) override;

  private:

    //! Donor material permittivity
    double _er;

    //! Donor exciton effective mass
    double _m;

    //! Forster radius
    double _Rf;

    //! Dexter radius
    double _Rd;

    //! Average donor-acceptor distance
    double _R_da;

    //! Total donor exciton lifetime
    double _tau;

    //! Donor exciton radiative lifetime
    double _tau_rad;

    //! Forster rate
    double _Kf;

    //! Dexter rate
    double _Kd;

    // donor and acceptor ids
    ID _id_d;
    ID _id_a;


};



//
// inline methods
// 

inline
ExcitonEnergyTransfer::ExcitonEnergyTransfer(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _Rf(3e-7),
    _Rd(1e-7),
    _R_da(1.5e-7),
    _er(1.0),
    _m(2.0),
    _tau(1e-9),
    _tau_rad(1e-9)
{
}





#endif // TC_EXCITONENERGYTRANSFER_H
