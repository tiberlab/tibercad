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
 * \file Alloy.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_ALLOY_H
#define TC_ALLOY_H


#include "tibercad/physics/Material.h"

//! Description of  an  alloy  material.
/*!
 * An alloy is thought to be of the form \f$A_xB_{1-x}C\f$.
 *
 * This class contains the properties of the constituent Materials A and B.
 */
class Alloy : public Material
{

  public:

    //! Destructor
    /*!
     * Deletes all \c PhysicalProperties objects
     */
    virtual ~Alloy(void);


    //! Create a material with name \c name
    static Alloy* create(const std::string& name,
        const ModelOptions& options);

    //! Get the molar fraction of the component A
    double get_molar_fraction(void) const;

    //! Get molar fraction for an atomic specie
    double get_molar_fraction(unsigned int label, Specie sp) const;

    //! Return the component A
    Material* get_component_A(void) const;

    //! Return the component B
    Material* get_component_B(void) const;


    //! Return the name of component material A
    const std::string& get_name_A(void) const;

    //! Return the name of component material B
    const std::string& get_name_B(void) const;

    //! Get the parent material for the pair of species (bond)
    const Material* get_parent(std::pair<Specie, Specie> atom_pair) const;

    //! decide whether an atom with a label can be substituted
    bool is_mutable(unsigned int) const; 

    //! get the map between species and molar fractions for a given label
    const std::map<Specie, double>& get_species_map(unsigned int) const;

    //! get the map between species and molar fractions for a given label
    const std::vector<std::map<Specie, double>>& get_species_map(void) const;

  protected:

    //! Construct an  alloy material
    Alloy(const std::string& name, const ModelOptions& options);


    //! \copydoc Material::do_preinit()
    virtual void do_preinit(void);


    //! \copydoc Material::do_init()
    virtual void do_init(void);


    //! Print info on Alloy composition
    virtual void do_info(void) const;


    //! set all species in the Alloy 
    virtual void fill_species(void);


  private:

    //! Molar fraction x of the \c Alloy
    /*!
     *  This is the \f$x\f$ in \f$A_xB_{1-x}C\f$
     */
    double _molar_fraction;


    //! The component A
    Material* _mat_A;


    //! The component B
    Material* _mat_B;

    //! Map for each label (starting from 1)
    /*! mapping Specie to molar fraction. E.g.
     *  In(x)Ga(1-x)As
     *  _specie_fraction[1]=<In,x> <Ga,1-x>
     *  _specie_fraction[2]=<As,1.0>
     */  
    std::vector<std::map<Specie, double>> _specie_fraction;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------



inline
const std::string&
Alloy::get_name_A(void) const
{
  return _mat_A->get_name();
}


inline
const std::string&
Alloy::get_name_B(void) const
{
  return _mat_B->get_name();
}


inline
double
Alloy::get_molar_fraction(void) const
{
  return _molar_fraction;
}

inline
double
Alloy::get_molar_fraction(unsigned int label, Specie sp) const
{
  std::map<Specie, double>::const_iterator it = _specie_fraction[label].find(sp);
  if (it == _specie_fraction[label].end() ) 
    return 0.0;
  else
    return it->second;
}


inline
Material*
Alloy::get_component_A(void) const
{
  return _mat_A;
}


inline
Material*
Alloy::get_component_B(void) const
{
  return _mat_B;
}

inline
const std::map<Specie, double>& 
Alloy::get_species_map(unsigned int i) const
{
  return _specie_fraction[i];
}

inline
const std::vector<std::map<Specie, double>>& 
Alloy::get_species_map(void) const
{
  return _specie_fraction;
}


#endif /* _ALLOY_H_ */
