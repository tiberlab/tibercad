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
 * \file Elementary1Diode.C
 * \brief tiberCAD pvmodule module implementation.
 *
 * \note This file is part of module pvmodule.
 */

#include "Elementary1Diode.h"
#include "Photocurrent.h"
#include "DegradationModel.h"

#include "tibercad/module/TiberModule.h"


Elementary1Diode::Elementary1Diode(const ModelOptions& options) :
  ElementaryCell(options)
{
}



void
Elementary1Diode::prepare_submodels(void)
{
  create_submodel(_photocurr_model, "photocurrent");
  create_submodel(_degradation_model, "degradation");
}

void
Elementary1Diode::do_init(void)
{
  _rseries = get_option("series_resistance", _rseries);
  _rshunt = get_option("shunt_resistance", _rshunt);
  _photocurr = get_option("photocurrent", _photocurr);
  _isat = get_option("saturation_current", _isat);
  _ideality = get_option("ideality_factor", _ideality);
}


void
Elementary1Diode::do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                  unsigned int& next_free,
                                  double area,
                                  const libMesh::Elem* elem,
                                  const libMesh::Point& p,
                                  std::ostream& os) const 
{
  double photocurr = _photocurr;
  double rseries = _rseries;
  double rshunt = _rshunt;
  double isat = _isat;
  double ideality = _ideality;

  if (_photocurr_model != nullptr)
    photocurr = _photocurr_model->get_photocurrent(elem, p);

  if (_degradation_model != nullptr)
  {
    DegradationModel::Parameters param;

    param.double_params.resize(5);
    // for now the DegradationModel and elementary cell model need to agree
    // implicitly on the order and meaning of parameters
    param.double_params[0] = rseries;
    param.double_params[1] = rshunt;
    param.double_params[2] = photocurr;
    param.double_params[3] = isat;
    param.double_params[4] = ideality;

    _degradation_model->degrade_parameters(elem, p, param);

    rseries = param.double_params[0];
    rshunt = param.double_params[1];
    photocurr = param.double_params[2];
    isat = param.double_params[3];
    ideality = param.double_params[4];
  }

  photocurr *= area;

  os << "Rs" << top_node << " " << top_node << " " << next_free << " " << rseries/area << "\n";
  os << "I" << top_node << " " << next_free << " " << bottom_node << " DC " << -photocurr << "\n";
  os << "D" << top_node << " " << next_free << " " << bottom_node << " DModel" << top_node << "\n";
  os << "Rsh" << top_node << " " << next_free << " " << bottom_node << " " << rshunt/area << "\n";
  os << ".model DModel" << top_node << " D(IS=" << isat*area << " N=" << ideality << ")\n";

  ++next_free;
}
