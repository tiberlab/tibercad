/*  
 * This file is part of the tiberCAD module pvmodule.
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
 * \file DegradationModel.h
 * \brief tiberCAD pvmodule module header.
 *
 * \note This file is part of module pvmodule.
 */

#ifndef TC_DEGRADATIONMODEL_H
#define TC_DEGRADATIONMODEL_H

#include "tibercad/physics/PhysicalModel.h"

/*!
 * \brief Base class for degradation models
 *
 * The model is intended for use with solar cell
 * equivalent circuits. It has to implement a method
 * that returns the model parameters changed according
 * to some degradation model.
 */
class DegradationModel : public PhysicalModel
{

  public:

    //! A base class for parameter containers
    /*!
     * Probably double values is what is needed, but
     * derived classes might need other stuff.
     */
    class Parameters {
      public:
        std::vector<double> double_params;
    };

    virtual ~DegradationModel(void) = default;

    //! Implementation for the 1-diode equivalent circuit
    void degrade_parameters(const libMesh::Elem* elem,
                            const libMesh::Point& p,
                            Parameters& params) const
    {
      do_degrade_params(elem, p, params);
    };


  protected:

    //! Constructor
    DegradationModel(const ModelOptions& options)
      : PhysicalModel(options) {};


    virtual void do_degrade_params(const libMesh::Elem* elem,
                                   const libMesh::Point& p,
                                   Parameters& params) const = 0;

  private:

};


#endif // TC_DEGRADATIONMODEL_H
