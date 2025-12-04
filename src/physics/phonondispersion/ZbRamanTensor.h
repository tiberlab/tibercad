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
 * \file ZbRamanTensor.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */

#ifndef _ZB_RAMANTENSOR_H_
#define _ZB_RAMANTENSOR_H_


#include "PhononModel.h"
#include "RamanTensor.h"
#include "Macrostrain.h"


class ZbRamanTensor: public RamanTensor
{
 public:
  //!constructor
  ZbRamanTensor(const ModelOptions& options) ;

  //!destructor
  ~ZbRamanTensor() {};

  //! Create a ZbLatticeThermalConductivity object
  static    ZbRamanTensor* create(const ModelOptions& options);

 
  virtual void re_init(){}; 

  virtual void set_phonon_model(PhononModel* phonon_model);
 
 private:

 PhononModel* _phonon_model;

 double _raman_d;
 
  protected:

  virtual void read_database(void);

  virtual void do_init(void);

  inline  virtual PhysicalModel*  create_new (void) const;


};

inline
ZbRamanTensor* ZbRamanTensor::create(const ModelOptions& options)
{
  return (new ZbRamanTensor(options));
}

inline
PhysicalModel*  ZbRamanTensor::create_new (void) const
{
  return (new  ZbRamanTensor(get_options()) );
}

inline
void  ZbRamanTensor::set_phonon_model(PhononModel* phonon_model)
{
  _phonon_model = phonon_model;
}



#endif
