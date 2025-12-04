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
 * \file WzDDsemiconductor.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef _WZDDSEMICONDUCTOR_H_
#define _WZDDSEMICONDUCTOR_H_


#include "tibercad/physics/semiconductormodels/DDsemiconductor.h"
#include "tibercad/physics/PhysicalModel.h"

#include<vector>

//! A class to provide all neccessary parameters for drift-diffusion calculation for a wurtzite  crystal.
class  WzDDsemiconductor : public DDsemiconductor
/*!
  The class can calculate information about the band structure, such as
  band edge energy, effective mass for the density of states calculation and
  degeneracy.
  Conduction band masses do not depend on strain.
  Only \f$ \Gamma \f$ minimum of conduction band is considered. 
*/
{

 public:

  

  //Destructor
  virtual ~WzDDsemiconductor(void) {};


  static WzDDsemiconductor* create(const ModelOptions& options);
 
 private:
 
  
  

 protected:

  //Constructor
  WzDDsemiconductor(const ModelOptions& options);

  PhysicalModel* create_new(void) const;
  
  //! calculates information about conduction bands
  /*!
    \f$ E_c^{\Gamma} = E_{c0}^{\Gamma} + a_{x}  (\varepsilon_{xx} + \varepsilon_{yy}) + a_z \varepsilon_{zz} \f$
  */
  virtual void  do_calculate_conduction_band_extremum(void);

};

inline
WzDDsemiconductor::WzDDsemiconductor(const ModelOptions& options)
 : DDsemiconductor(options)
{
}

inline PhysicalModel* WzDDsemiconductor::create_new( ) const
{
  return ( new WzDDsemiconductor(get_options()) );
}

inline WzDDsemiconductor* WzDDsemiconductor::create(const ModelOptions& options)
{
  return  new WzDDsemiconductor(options) ;
}





#endif 
