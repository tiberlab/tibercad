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
 * \file Dopant.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _DOPANT_H_
#define _DOPANT_H_

#include "tibercad/base/ModelOptions.h"
#include "tibercad/base/tiber_dll.h"

#include <string>

class ExternalProfile;

namespace libMesh
{
class Elem;
class Point;
}


//! Describes a dopant with a single energy level
class TBDLEXPORT Dopant
{

  public:


    //! The type of the dopant
    enum DopingType
    {
      P_TYPE = -1, /*!< p-type, acceptor */
      N_TYPE = 1   /*!< n-type, donor */
    };

    //! Constructor
    /*!
     * \param density doping density
     * \param ionisation_energy the distance of the doping level from
     * the corresponding band edge \f$\vert E_{c,v} - E_d\vert\f$
     * \param g_factor g (cf. get_ionized_dopant_density())
     * \param type the doping type
     */
    Dopant(double density = 0.0, double ionisation_energy = 0.025,
        int g_factor = 2, DopingType type = N_TYPE);


    //! Copy constructor
    Dopant(const Dopant& dopant);


    // Destructor
    virtual ~Dopant(void) {};


    //! Create a doping with given profile and options
    static Dopant* create(const std::string& profile, const ModelOptions& options);

    
    //! Get the doping density
    double get_doping_density(void) const;

    
    //! Get the ionisation energy of the dopant
    /*!
     * \f$E_i = \vert E_{c,v} - E_d\vert\f$
     */
    double get_ionisation_energy(void) const;

    
    //! Get g (cf. get_ionized_dopant_density())
    int get_g_factor(void) const;

    
    //! Get the doping type
    DopingType get_type(void) const;

    
    //! Calculate the coordinate-dependent doping density
    /*! 
     * Sets the local doping density
     */
    void calculate_doping_density(const libMesh::Elem* elem, const libMesh::Point& p);


    //! Get the density of ionized dopants
    /*!
     * For donors this is
     * \f[
     * N_d^+ = \frac{N_d}{1 + g \exp{??}}
     * \f]
     *
     */
    double get_ionized_dopant_density(double arg, double kT);


    /*! \brief {Get the derivative of the ionized doping density with
     * respect to the electric potential}
     */
    double get_ionized_dopant_density_derivative(double arg, double kT);

    //! Get the options
    ModelOptions& get_options(void);



  protected:

    //! Calculate the doping density
    virtual double do_calculate_doping_density(const libMesh::Elem* elem, const libMesh::Point& p);



  private:

    //! The doping profile function
    ExternalProfile* _profile;

    //! The doping density
    double _density;

    //! The doping type
    DopingType _type;

    //! Whether to use incomplete ionisation
    bool _incomplete_ionisation;

    //! The ionisation energy
    double _ionisation_energy;

    //! The g factor
    int _g_factor;


    //! More options
    ModelOptions _options;

};



//
// inline member functions
//

inline
Dopant::Dopant(double density, double ionisation_energy,
               int g_factor, DopingType type)
  : _profile(nullptr),
    _density(density),
    _type(type),
    _incomplete_ionisation(true),
    _ionisation_energy(ionisation_energy),
    _g_factor(g_factor)
{
}

inline
Dopant::Dopant(const Dopant& dopant)
  : _profile(dopant._profile),
    _density(dopant._density),
    _type(dopant._type),
    _incomplete_ionisation(dopant._incomplete_ionisation),
    _ionisation_energy(dopant._ionisation_energy),
    _g_factor(dopant._g_factor)
{
}

inline
double
Dopant::get_doping_density(void) const
{
  return _density;
}



inline
double
Dopant::get_ionisation_energy(void) const
{
  return _ionisation_energy;
}

inline
int
Dopant::get_g_factor(void) const
{
  return _g_factor;
}

inline
Dopant::DopingType
Dopant::get_type(void) const
{
  return _type;
}


inline
ModelOptions&
Dopant::get_options(void)
{
  return _options;
}



inline
double
Dopant::do_calculate_doping_density(const libMesh::Elem*, const libMesh::Point&)
{
  return 0.0;
}



#endif //_DOPANT_H_
