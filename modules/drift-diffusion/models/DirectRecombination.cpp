/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file DirectRecombination.C
 * \brief tiberCAD driftdiffusion module implementation.
 *
 * \note This file is part of module driftdiffusion.
 */


#include "DirectRecombination.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/io/Database.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/io/Messages.h"

#include "mesh_base.h"
#include "quadrature.h"
#include "enum_quadrature_type.h"

#include "tibercad/module/TiberModule.h"


using namespace std;

DirectRecombination::QRecMap
DirectRecombination::_qrec_vals;


DirectRecombination::DirectRecombination(const ModelOptions& options)
  : RecombinationModelInterface(options),
    C_(0.0),
    _quantum_optics(NULL)
{
  set_name("RadiativeRecombination");
}


void
DirectRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/direct");

  C_ = db.get("C", C_);
  _gamma = db.get("gamma", _gamma);

}



void
DirectRecombination::do_init(void)
{
  get_parameter("C", C_);
  get_parameter("gamma", _gamma);

  string quantumsim = get_option("optics_simulation", "");
  if (!quantumsim.empty())
  {
    _quantum_optics = SimulationInterface::find_simulation(quantumsim);
    if (_quantum_optics == NULL)
      throw InitFailedException("Cannot find optics simulation \'" + quantumsim + "\'");

    _rec_id = _quantum_optics->get_solution_id("Recombination");
    if (_rec_id == INVALID_ID)
      throw InitFailedException("Simulation \'" + quantumsim + "\'" +
          " does not have the needed solution \'Recombination\'");

    // TODO should check for consistency of regions
  }

}



void
DirectRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double T = dd.get_lattice_temperature();

  double C = C_*pow(T/0.02585, _gamma);

  double f = 1.0 - exp((Efp - Efn) / T);
  recomb_e = recomb_h = C * n * p * f;
}



void
DirectRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double Efn  = -dd.get_electron_electro_chemical_potential();
  double Efp  = -dd.get_hole_electro_chemical_potential();
  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double T = dd.get_lattice_temperature();

  double C = C_*pow(T/0.02585, _gamma);

  double expf = exp((Efp - Efn) / T);
  double f = 1.0 - expf;

  recomb_e[0] = recomb_h[0] = C * p * f; // dR/dn
  recomb_e[1] = recomb_h[1] = C * n * f; // dR/dp
  recomb_e[2] = recomb_h[2] = -C * n * p / T * expf; // dR/dEfn
  recomb_e[3] = recomb_h[3] = C * n * p / T * expf; // dR/dEfn
}


void
DirectRecombination::do_reinit(void)
{
  if (_quantum_optics != NULL)
  {
    SimulationInterface* sim = SimulationInterface::get_simulation(get_simulator_id());
    unsigned int seqnum = sim->get_solve_sequence_number();

    QRecMap::iterator it(_qrec_vals.find(make_pair(_quantum_optics, sim)));
    if ((it != _qrec_vals.end()) && ((it->second).first == seqnum))
      C_ = (it->second).second;
    else
    {
      map<ID, vector<double> > data;
      data[_rec_id];

      if (_quantum_optics->get_solution(data))
      {
        double rec = data[_rec_id][0];
        Messages::info("Recalculate radiative recombination parameter: ", false);

        // now we have to integrate the term (np - ni^2)
        // for that, we have to loop over all elements
        // TODO this has to be checked for parallel execution

        data.clear();
        // TODO IntrinsicDensity is currently missing
        ID edens_id = sim->get_solution_id("eDensity");
        ID hdens_id = sim->get_solution_id("hDensity");
        //ID idens_id = sim->get_solution_id("IntrinsicDensity");
        data[edens_id];
        data[hdens_id];

        unsigned int dim = sim->get_mesh().mesh_dimension();
        std::unique_ptr<libMesh::FEBase> fe(sim->build_finite_element(dim, libMesh::FEType()));
        std::unique_ptr<libMesh::QBase> qrule(libMesh::QBase::build(libMesh::QGAUSS, dim, libMeshEnums::FIFTH));
        fe->attach_quadrature_rule(qrule.get());

        const vector<Real>& JxW = fe->get_JxW();

        double tot_rec = 0.0;

        MeshBase::const_element_iterator it = sim->active_local_elements_begin();
        const MeshBase::const_element_iterator end = sim->active_local_elements_end();
        for ( ; it != end; ++it)
        {
          const libMesh::Elem* elem = *it;
          if (_quantum_optics->includes_region(elem->subdomain_id()))
          {
            fe->reinit(elem);

            if (sim->get_solution(elem, data, qrule->get_points(), true))
            {
              for (int n = 0; n < qrule->n_points(); ++n)
              {
                double np = data[edens_id][n] * data[hdens_id][n];
                tot_rec += JxW[n] * np;
              }
            }
          }
        }

        C_ = rec / tot_rec;
        _qrec_vals[make_pair(_quantum_optics, sim)] = make_pair(seqnum, C_);

        ostringstream os;
        os << "B = " << rec / tot_rec << endl;
        Messages::info(os.str());
        Messages::newline();
      }
    }
  }
}
