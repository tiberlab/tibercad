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
 * \file ZbFreeDynamicalMatrix.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */

#ifndef TC_ZB_FREEDYNAMICALMATRIX_H
#define TC_ZB_FREEDYNAMICALMATRIX_H


#include "PhononModel.h"
#include "DynamicalMatrix.h"

class ZbFreeDynamicalMatrix: public DynamicalMatrix
{
 public:
  //!constructor
  ZbFreeDynamicalMatrix(const ModelOptions& options) : DynamicalMatrix(options) {};

  //!destructor
  ~ZbFreeDynamicalMatrix() {};

  //! Create a ZbLatticeThermalConductivity object
 static  ZbFreeDynamicalMatrix* create(const ModelOptions& options);


//! Update the lattice thermal conductivity given the Temperature
  virtual void re_init(); 


  virtual void set_phonon_model(PhononModel* phonon_model);
 
 private:

 PhononModel* _phonon_model;
 double w0;

 protected:

  virtual void read_database(void);

  virtual void do_init(void);

  inline  virtual PhysicalModel*  create_new (void) const;

};

inline
ZbFreeDynamicalMatrix* ZbFreeDynamicalMatrix::create(const ModelOptions& options)
{
  return (new ZbFreeDynamicalMatrix(options));
}

inline
PhysicalModel*  ZbFreeDynamicalMatrix::create_new (void) const
{
  return (new  ZbFreeDynamicalMatrix(get_options()) );
}

inline
void ZbFreeDynamicalMatrix::set_phonon_model(PhononModel* phonon_model)
{
  _phonon_model = phonon_model;
}



#endif
