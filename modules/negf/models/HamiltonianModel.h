/*  
 * This file is part of the tiberCAD module negf.
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
 * \file HamiltonianModel.h
 * \brief tiberCAD negf module header.
 *
 * \note This file is part of module negf.
 */


#include "tibercad/physics/PhysicalModel.h"
#include "tensor_value.h"
#include "vector_value.h"

class HamiltonianModel : public PhysicalModel
{

  public:

    virtual ~HamiltonianModel(void) {};

    const libMesh::TensorValue<double>& get_inv_mass(void) const;

    double get_degeneracy(void) const;

    const std::string& get_band_type(void) const;

    const std::string& get_model_name(void) const;
    
    const std::string& get_simulation(void) const;

    int get_n_bands(void) const;

  protected:

    HamiltonianModel(const ModelOptions& options);

    libMesh::TensorValue<double> _inv_mass_crys;

    libMesh::TensorValue<double> _inv_mass;

    double _degeneracy;

    int _num_bands;

    std::string _band_type;

    std::string _model;

    std::string _simulation;
};


inline
const libMesh::TensorValue<double>&
HamiltonianModel::get_inv_mass(void) const
{
  return _inv_mass;
}

inline
double
HamiltonianModel::get_degeneracy(void) const
{
  return _degeneracy;
}

inline
const std::string&
HamiltonianModel::get_band_type(void) const
{
  return _band_type;
}

inline
const std::string&
HamiltonianModel::get_model_name(void) const
{
  return _model;
}

inline
const std::string&
HamiltonianModel::get_simulation(void) const
{
  return _simulation;
}

inline
int
HamiltonianModel::get_n_bands(void) const
{
  return _num_bands;
}

inline
HamiltonianModel::HamiltonianModel(const ModelOptions& options)
  : PhysicalModel(options)
{
}


