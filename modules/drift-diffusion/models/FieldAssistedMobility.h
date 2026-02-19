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
 * \file FieldAssistedMobility.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_FIELDASSISTEDMOBILITY_H
#define TC_FIELDASSISTEDMOBILITY_H

#include "MobilityModelInterface.h"


//! Field-assisted mobility model for unordered systems
/*!
 * The mobility is calculated as
 * \f[
 * \mu = \mu_0 e^{\sqrt{|E|/E_0}}
 * \f]
 */
class TC_DLLOCAL FieldAssistedMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~FieldAssistedMobility(void) {};

    //! Create a FieldAssistedMobility object
    static FieldAssistedMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm);


  protected:

    //! constructor
    FieldAssistedMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModel* create_new(void) const;


  private:

    enum DrivingForce
    {
      EFIELD,
      GRADFERMI,
    };


    //! The zero-field mobility
    double _mu0;


    //! The critical field strength
    double _E0;


    //! The driving force to be used
    DrivingForce _force;


};

//
// inline methods
//

inline
FieldAssistedMobility::FieldAssistedMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    _mu0(0.0054),
    _E0(3e5),
    _force(EFIELD)
{
}


inline
FieldAssistedMobility*
FieldAssistedMobility::create(const ModelOptions& options)
{
  return new FieldAssistedMobility(options);
}


inline
PhysicalModel*
FieldAssistedMobility::create_new(void) const
{
  return new FieldAssistedMobility(get_options());
}


#endif // TC_FIELDASSISTEDMOBILITY_H
