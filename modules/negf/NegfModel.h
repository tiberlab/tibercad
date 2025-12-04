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
 * \file NegfModel.h
 * \brief tiberCAD negf module header.
 *
 * \note This file is part of module negf.
 */

/*
 * Negfmodel.h
 *
 *  Created on: Jan 18, 2012
 *      Author: fpalomba
 */

#ifndef NEGFMODEL_H_
#define NEGFMODEL_H_

#include "tibercad/physics/PhysicalModel.h"
#include "models/HamiltonianModel.h"

class TBDLLOCAL NegfModel : public PhysicalModel
{
  public:

    //! Destructor
    virtual ~NegfModel(void) {};

    //! Creator function
    static NegfModel* create(const Material* mat, const ModelOptions& options);

    const libMesh::TensorValue<double>& get_inv_mass(ID id) const;

    unsigned int get_n_bands(void) const;

    const std::string& get_band(ID id) const;

    double get_degeneracy(ID id) const;

    std::string get_simulation(ID id) const;
    
    std::string get_model_name(ID id) const;

    HamiltonianModel* get_Hamiltonian_model(ID id) const;

  protected:

    virtual void prepare_submodels(void);

  private:

    static TiberModelObject*  _create(const ModelOptions& options, const void*);

    static void  _destroy( TiberModelObject* p);

    NegfModel(const ModelOptions& options);

    std::vector<HamiltonianModel*> _ham_models;
};

inline
TiberModelObject*  NegfModel::_create(const ModelOptions& options, const void*)
{
  return new NegfModel(options);
}

inline
void  NegfModel::_destroy( TiberModelObject* p)
{
  delete p;
}


inline
const libMesh::TensorValue<double>& NegfModel::get_inv_mass(ID id) const
{
   return _ham_models[id]->get_inv_mass();
}

inline
unsigned int NegfModel::get_n_bands(void) const
{
  int n_bands = 0;

  for (unsigned int i=0; i<_ham_models.size(); i++)
  {
    n_bands += _ham_models[i]->get_n_bands();
  }

  return n_bands;
}

inline
const std::string& NegfModel::get_band(ID id) const
{
   return _ham_models[id]->get_band_type();
}

inline
double NegfModel::get_degeneracy(ID id) const
{
   return _ham_models[id]->get_degeneracy();
}

inline
std::string NegfModel::get_simulation(ID id) const
{
   return _ham_models[id]->get_simulation();
}

inline
std::string NegfModel::get_model_name(ID id) const
{
   return _ham_models[id]->get_model_name();
}

inline
HamiltonianModel* NegfModel::get_Hamiltonian_model(ID id) const
{
  return _ham_models[id];
}

#endif /* NEGFMODEL_H_ */
