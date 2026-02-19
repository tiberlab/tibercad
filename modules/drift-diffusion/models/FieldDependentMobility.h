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
 * \file FieldDependentMobility.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_FIELDDEPENDENTMOBILITY_H
#define TC_FIELDDEPENDENTMOBILITY_H

#include "MobilityModelInterface.h"


//! FieldDependent mobility model
/*!
 * The constant mobility model assumes no mobility dependence on
 * doping density or electric field. Only temperature dependence is
 * included.
 * The mobility is calculated from
 * \f[
 * \mu = \mu_{max} \left(\frac{T}{T_0}\right)^{-\gamma}
 * \f]
 */
class TC_DLLOCAL FieldDependentMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~FieldDependentMobility(void);

    //! Create a FieldDependentMobility object
    static FieldDependentMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm);


  protected:

    //! constructor
    FieldDependentMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! Create low field mobility model
    virtual void prepare_submodels(void);



  private:

    enum DrivingForce
    {
      EFIELD,
      GRADFERMI,
      FIELDPARAM
    };


    //! The exponent
    double _beta;


    //! The exponent for the temperetaure dependence of _beta
    double _betaexp;


    //! The maximum saturation velocity
    double _vsat0;


    //! For the temparature dependence of v_sat
    double _vsat_b;


    //! The minimum of vsat for formula 2
    double _vsat_min;


    //! The formula to be used
    int _vsat_formula;


    //! The low-field mobility
    MobilityModelInterface* _low_field_mob;


    //! The driving force to be used
    DrivingForce _force;


    //! A damping parameter
    double _damping;


};

//
// inline methods
//

inline
FieldDependentMobility::FieldDependentMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    _beta(1),
    _betaexp(0.0),
    _vsat0(1.13e7),
    _vsat_b(1),
    _vsat_min(5e5),
    _vsat_formula(1),
    _low_field_mob(NULL),
    _force(EFIELD),
    _damping(1e9)
{
}


inline
FieldDependentMobility*
FieldDependentMobility::create(const ModelOptions& options)
{
  return new FieldDependentMobility(options);
}



inline
FieldDependentMobility::~FieldDependentMobility(void)
{
}

#endif // TC_FIELDDEPENDENTMOBILITY_H
