/*  
 * This file is part of the tiberCAD module masstransport.
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
 * \file MTModel.h
 * \brief tiberCAD masstransport module header.
 *
 * \note This file is part of module masstransport.
 */


#ifndef TC_POISSONMODEL_H
#define TC_POISSONMODEL_H

#include "tibercad/physics/PhysicalModel.h"


//! This is the base class for the MT physical model
class TC_DLEXPORT MTModel : public PhysicalModel
{

  public:

    //! Destructor
    virtual ~MTModel(void);

    //! Creator function
    static MTModel* create(const Material* mat, const ModelOptions& options);

    //! Calculate everything
    void calculate(const Elem* elem, const Point& point);

    //! Get the solubility in g/m^3/P
    double get_solubility(void) const;

    //! Get the diffusivity in m^2/s
    double get_diffusivity(void) const;


  protected:

    //! Constructor
    MTModel(const ModelOptions& options);

    virtual void do_init(void) override;

    virtual void prepare_submodels(void) override;

  

  private:
 
    //! The solubility in g/m^3/P
    double _solubility = 0.45; // this value is for H2O
 
    //! The diffusivity in m^2/s
    double _diffusivity = 3.43e-11; // this value is for H2O

    //! The constructor method
    static TiberModelObject* _create(const ModelOptions& options);

    //! The destructor method
    static void _destroy(TiberModelObject* p);

    //! The cell temperature
    double _cell_temp = 300;

};




inline
MTModel::MTModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


inline
double
MTModel::get_solubility(void) const
{
  return _solubility;
}

inline
double
MTModel::get_diffusivity(void) const
{
  return _diffusivity;
}


#endif // TC_POISSONMODEL_H



