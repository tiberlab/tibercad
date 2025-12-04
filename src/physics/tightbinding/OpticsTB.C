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
 * \file OpticsTB.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/tightbinding/OpticsTB.h"
#include "tibercad/physics/schroedinger/EigenvalueProblem.h"
#include "tibercad/io/DataOutput.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/io/Messages.h"

#include "mesh.h"
#include "elem.h"
#include "mesh_generation.h"

using namespace Constants;


OpticsTB::OpticsTB(const ModelOptions& options)
 : Optics(options)
{
}


OpticsTB::~OpticsTB()
{
}


//===============================================//
void OpticsTB::do_init()
{

  std::cout<<"tb optics init"<<std::endl;

  Optics::do_init();

  std::cout<<"tb optics initialized"<<std::endl;

}

//===============================================//




//=========================================================================================

void OpticsTB::do_compute_matrix_elements( )
{
  unsigned int n_i =  _initial_indices.size();
  unsigned int n_f =  _final_indices.size();

  for (unsigned int i = 0; i < 3; i++)
  {
    _P_matrix[i].resize(n_i);
    for (unsigned int j = 0; j < n_i; j++)  _P_matrix[i][j].resize(n_f);
  }

  ModelOptions options;

  for (unsigned int i = 0; i < 3; i++)
  {
     options.set_option("P_matrix", true);
     options.set_option("poldir", i+1);

     assemble(options);


     for (unsigned int i1 = 0; i1 < n_i; i1++)
       for (unsigned int i2 = 0; i2 < n_f; i2++)
       {

          unsigned int is = _initial_state_numbers[_initial_indices[i1]];
          unsigned int fs = _final_state_numbers[_final_indices[i2]];

          _P_matrix[i][i1][i2] = 
              _initial_state_model->calculate_matrix_element(
                             _initial_state_particle, is,
		                         _final_state_particle, fs);

          // convert matrix element from eV*Ang to atomic units
          // 1 eV = 1 H / 27.2114
          // 1 Ang = 1/0.529177 a0 
          // eV*Ang = 1 H*a0 /(27.2114 * 0.529177) 
          // 1 bhor_rad = 5.2917721e-11 m * 10^10 Ang/m
          _P_matrix[i][i1][i2] *= 1/(Hartree*bohr_radius*1e10);

       }
   }


}

//=========================================================================================
void
OpticsTB::do_assemble(const ModelOptions& options)
{
  _initial_state_model->assemble(options);
}	  

