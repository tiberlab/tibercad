/*  
 * This file is part of the tiberCAD module wateringress.
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
 * \file WIModel.h
 * \brief tiberCAD wateringress module header.
 *
 * \note This file is part of module wateringress.
 */


#ifndef TC_POISSONMODEL_H
#define TC_POISSONMODEL_H

#include "tibercad/physics/PhysicalModel.h"


//! This is the base class for the WI physical model
class TC_DLEXPORT WIModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~WIModel(void);

    //! Creator function
    static WIModel* create(const Material* mat, const ModelOptions& options);

    //! Calculate everything
    void calculate(const Elem* elem, const Point& point);

    //! Get the water solubility in g/m^3/P
    double get_solubility(void) const;

    //! Get the water diffusivity in m^2/s
    double get_diffusivity(void) const;


  protected:

    //! Constructor
    WIModel(const ModelOptions& options);

    virtual void do_init(void) override;

    virtual void prepare_submodels(void) override;

  

  private:
 
    //! The water solubility in g/m^3/P
    double _solubility = 0.45;
 
    //! The water diffusivity in m^2/s
    double _diffusivity = 3.43e-11;

    //! The constructor method
    static TiberModelObject* _create(const ModelOptions& options);

    //! The destructor method
    static void _destroy(TiberModelObject* p);

};




inline
WIModel::WIModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


inline
double
WIModel::get_solubility(void) const
{
  return _solubility;
}

inline
double
WIModel::get_diffusivity(void) const
{
  return _diffusivity;
}


#endif // TC_POISSONMODEL_H
