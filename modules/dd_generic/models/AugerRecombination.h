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
 * \file AugerRecombination.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef _AUGERRECOMBINATION_H_
#define _AUGERRECOMBINATION_H_

#include "RecombinationModelInterface.h"
#include "tibercad/base/TypeDefs.h"


//! Implementation of Auger recombination
/*!
 * This class implements Auger recombination processes that can be
 * modeled by 
 * \f[
 *   R_{aug} = (C_nn + C_pp)(np-n_i^2)
 * \f]
 * with
 * \f[
 *  C_{\{n,p\}} = \left(A + B\frac{T}{T_0} + C\left(\frac{T}{T_0}\right)^2\right)
 *      \left(1 + H e^{-\{n,p\}/N_0}\right)
 * \f]
 */
class TBDLLOCAL AugerRecombination : public RecombinationModelInterface
{

  public:

    //! Destructor
    virtual ~AugerRecombination(void) {};

    //! Create a ConstantMobility object
    static AugerRecombination* create(const ModelOptions& options);

    
  protected:

    //! Constructor
    AugerRecombination(const ModelOptions& options);

    //! \copydoc RecombinationModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc RecombinationModelInterface::do_init()
    virtual void do_init(void) override;

    //! \copydoc RecombinationModelInterface::calculate_rate_and_derivatives()
    virtual void calculate_rate_and_derivatives(std::vector<double>& R,
        std::vector<std::vector<double>>& dPotentials) override;

    
  private:

    //! The parameter A for the electrons
    double _An;

    //! The parameter A for the holes
    double _Ap;

    //! The parameter B for the electrons
    double _Bn;

    //! The parameter B for the holes
    double _Bp;

    //! The parameter C for the electrons
    double _Cn;

    //! The parameter C for the holes
    double _Cp;

    //! The parameter H for the electrons
    double _Hn;

    //! The parameter H for the holes
    double _Hp;

    //! The parameter N0 for the electrons
    double _N0n;

    //! The parameter N0 for the holes
    double _N0p;


    //! Whether to use fixed Cn or not
    bool _fixed_Cn;

    //! Whether to use fixed Cp or not
    bool _fixed_Cp;


    //! Get Cn
    double get_Cn(void);

    //! Get Cp
    double get_Cp(void);

};


//
// inline methods
// 


inline
AugerRecombination::AugerRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _An(6.7000e-32),
    _Ap(7.2000e-32),
    _Bn(2.4500e-31),
    _Bp(4.5000e-33),
    _Cn(-2.2000e-32),
    _Cp(2.6300e-32),
    _Hn(3.46667),
    _Hp(8.25688),
    _N0n(1.0000e+18),
    _N0p(1.0000e+18),
    _fixed_Cn(false),
    _fixed_Cp(false)
{
}


inline
AugerRecombination*
AugerRecombination::create(const ModelOptions& options)
{
  return new AugerRecombination(options);
}





#endif // _AUGERRECOMBINATION_H_
