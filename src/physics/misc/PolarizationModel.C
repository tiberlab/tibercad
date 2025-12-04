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
 * \file PolarizationModel.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/misc/PolarizationModel.h"
#include "tibercad/atomistic/BulkCrystal.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/base/InitFailedException.h"
#include "tibercad/io/Messages.h"


using namespace std;


PolarizationModel::PolarizationModel(const ModelOptions& options) :
  PhysicalModel(options),
  _polarization(0),
  _polarization_source(NULL, INVALID_ID),
  _fixed_or_external(false)
{
}


PolarizationModel*
PolarizationModel::create(const ModelOptions& options)
{
  return(new PolarizationModel(options));
}


void
PolarizationModel::calculate(const libMesh::Elem* elem, const libMesh::Point& point)
{
  if (!_fixed_or_external)
    do_calculate(elem, point);
  else
    PolarizationModel::do_calculate(elem, point);
}


void 
PolarizationModel::rotate(void)
{

  if (!_fixed_or_external)
  {
    const Material* mat = get_material();
    _polarization = mat->get_rotation_matrix() * _polarization;
  }

}


void
PolarizationModel::do_print_info(void)
{
  if (_fixed_or_external)
  {
    Messages::info("Simple polarization model");
  }
}




void
PolarizationModel::do_calculate(const Elem* elem, const Point& point)
{
  if (_polarization_source.first != NULL)
  {
    vector<double> values(3);
    vector<Point> p(1, point);
    ID id = _polarization_source.second;

    if (_polarization_source.first->get_solution(elem, id, values, p))
    {
      _polarization(0) = values[0];
      _polarization(1) = values[1];
      _polarization(2) = values[2];
    }
  }
}


void
PolarizationModel::do_init(void)
{
  string polarization_str = get_option("polarization", "");
  vector<double> polarization;
  Utils::extract_vector(polarization_str, polarization);

  if (polarization.size() == 3)
  {
    // seems to be a polarization vector
    _polarization(0) = polarization[0];
    _polarization(1) = polarization[1];
    _polarization(2) = polarization[2];
    _fixed_or_external = true;
  }
  else if (!polarization_str.empty())
  {
    _polarization_source =
      SimulationInterface::find_solution_provider(polarization_str);
    if ((_polarization_source.first == NULL) || (_polarization_source.second == INVALID_ID))
    {

      throw InitFailedException(polarization_str + " is invalid identifier for "
          "a module providing polarization field.");
    }
    _fixed_or_external = true;
  }
}
