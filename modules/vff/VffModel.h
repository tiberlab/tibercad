/*  
 * This file is part of the tiberCAD module vff.
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
 * \file VffModel.h
 * \brief tiberCAD vff module header.
 *
 * \note This file is part of module vff.
 */

#ifndef TC_VFFMODEL_H
#define TC_VFFMODEL_H


#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/io/Messages.h"
#include "tibercad/io/Database.h"
#include "Keating.h"
#include "tibercad/atomistic/Atom.h"

class Keating;


//! This is the base class for VFF models
/*!
 * Right now there are no derived classes from this, so it is declared local,
 * but in future there might be a need to specialize this class.
 */
class TC_DLLOCAL VffModel : public PhysicalModel
{

public:

  //! Destructor
  virtual ~VffModel(void);


  //! Creator function
  static VffModel* create(const Material* mat, const ModelOptions& options);

  double get_alpha(void) const;

  double get_beta(void) const;

  double get_costeta(void) const;

  double get_d(void) const;

  double get_alpha(const Atom& atm1, const Atom& atm2) const;

  double get_beta(const Atom& atm1, const Atom& atm2, const Atom& atm3) const;

  double get_costeta(const Atom& atm1, const Atom& atm2, const Atom& atm3) const;

  double get_d(const Atom& atm1, const Atom& atm2) const;


protected:

  //! Constructor
  VffModel(const ModelOptions& options);

  //! Read database
  void read_database(void);

  //! Initialize
  virtual void do_init(void);

  //virtual void prepare_submodels(void);
  virtual void do_print_info(void);

private:


  static TiberModelObject* _create(const ModelOptions& options);

  static void  _destroy( TiberModelObject* p);

  void prepare_submodels(void);

  Keating* _keating;

  bool along_c(const Atom& atm1, const Atom& atm2) const;

};


inline
VffModel::VffModel(const ModelOptions& options) :
PhysicalModel(options)
{
}

inline
VffModel::~VffModel()
{
}

inline
TiberModelObject*  VffModel::_create(const ModelOptions& options)
{

  return new VffModel(options);

}

inline
void  VffModel::_destroy( TiberModelObject* p)
{

  delete p;

}

inline
double VffModel::get_alpha(void) const
{
  return _keating->get_alpha_0();
}

inline
double VffModel::get_beta(void) const
{
  return _keating->get_beta_0();
}

inline
double VffModel::get_costeta(void) const
{
  return _keating->get_costeta_0();
}

inline
double VffModel::get_d(void) const
{
  return _keating->get_d_0();
}

#endif // TC_VFFMODEL_H
