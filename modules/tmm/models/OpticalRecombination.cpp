/*  
 * This file is part of the tiberCAD module tmm.
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
 * \file OpticalRecombination.C
 * \brief tiberCAD tmm module implementation.
 *
 * \note This file is part of module tmm.
 */


#include "OpticalRecombination.h"
#include "tibercad/module/TiberModule.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/io/Messages.h"
#include <string>



void 
OpticalRecombination::calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda)
{
  string recombination_str(get_option("radiative_recombination", "0"));
   vector<string> recombination;
   Utils::extract_vector(recombination_str, recombination);
   _recombination_model.resize(recombination.size());
   _recombination_id.resize(recombination.size());
   for (size_t i = 0; i < recombination.size(); ++i)
    {
      pair<SimulationInterface*, ID> provider =
        SimulationInterface::find_solution_provider(recombination[i]);
      _recombination_model[i] = provider.first;
      _recombination_id[i] = provider.second;
      if (_recombination_model[i] == NULL)
        throw InitFailedException("Cannot find radiative_recombination model: " + recombination[i]);
    }
	
    _recombination_rate = 0;
	

    vector<double> tmp(1);
    for (size_t i = 0; i < _recombination_model.size(); ++i)
    {
      if (_recombination_model[i]->get_solution(elem, _recombination_id[i], tmp,
          vector<Point>(1, point)))
        _recombination_rate += tmp[0];
    }

	

	set_emission_power(_recombination_rate ); 
}

