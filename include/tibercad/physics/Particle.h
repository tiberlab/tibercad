/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file Particle.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



//! A class to represent a particle population
/*!
 * The particle population is represented by the 
 * type of particle (its charge), its density,
 * an equivalent quasi-Fermi level according to 
 * a thermal equilibrium Fermi-Dirac statistics,
 * an effective carrier temperature and by the
 * derivative with respect to the quasi Fermi level.
 *
 * This class is mainly thought for the exchange of
 * density data between models.
 */
class Particle
{

  public:

    //! Constructor
    Particle(double charge, double density, double Ef, double kT);

    //! Get the density
    double density(void) const
    { return _density; }

    //! Get the density
    double kT(void) const
    { return _kT; }

    //! Get the quasi Fermi level
    double fermi_level(void) const
    { return _fermi_level; }


  private:

    //! The charge in units of the elementary charge
    double _charge;

    //! The density in cm^-3
    double _density;

    //! The quasi Fermi level
    double _fermi_level;

    //! The effective temperature in eV
    double _kT;

};

inline
Particle::Particle(double charge, double density, double Ef, double kT) :
  _charge(charge),
  _density(density),
  _fermi_level(Ef),
  _kT(kT)
{
}
