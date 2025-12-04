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
 * \file ZbDDsemiconductor.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _ZBDDSEMICONDUCTOR_H_
#define _ZBDDSEMICONDUCTOR_H_
 
 
#include "tibercad/physics/semiconductormodels/DDsemiconductor.h"
#include "tibercad/physics/PhysicalModel.h"
#include<vector>
#include<complex>

 
class ZbDDsemiconductor  : public DDsemiconductor
{
 public:
  
   
 
  
  virtual ~ZbDDsemiconductor(void) {};
 
  
 
  static ZbDDsemiconductor* create(const ModelOptions& options);
  
 private:
 
 
 
 protected:
 
  ZbDDsemiconductor(const ModelOptions& options) : DDsemiconductor(options) {};

  virtual PhysicalModel* create_new(void) const;
   
  
  virtual void  do_calculate_conduction_band_extremum(void);
   
  
};
 
 
 
inline PhysicalModel* ZbDDsemiconductor::create_new( ) const
{
  return ( new ZbDDsemiconductor(get_options()) );
}
 
inline ZbDDsemiconductor* ZbDDsemiconductor::create(const ModelOptions& options)
{
  return new ZbDDsemiconductor(options);
}
 
#endif


