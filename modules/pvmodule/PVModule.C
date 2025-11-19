// $Id: PVModule.C 5514 2024-05-03 15:22:47Z maufder $

#include "PVModule.h"
#include "PVModuleModel.h"
#include "PVModuleBoundaryModel.h"
#include "ElementaryCell.h"
#include "tibercad/solver/TiberLinearSystem.h"
#include "tibercad/io/Messages.h"
#include "tibercad/io/Database.h"

#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"

// This is needed in order to create the shared module library
#include "tibercad/module/TiberModule.h"

#include <fstream> 
#include <sstream>
#include <queue>

using namespace std;
using namespace libMesh;



PVModule::PVModule(const ModelOptions& options) :
  SimulationInterface(options),
  _my_assembly(this)
{
  // there's nothing to be done
}


PVModule::~PVModule(void)
{
  // there's nothing to be done
}


PVModule*
PVModule::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new PVModule(options);
}



void
PVModule::do_init(void)
{
  parse_options();

  // create a linear equation system 
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // add variables and attach the assemble function
  if (_discretization == DEC)
  {
    system.add_variable("Vtop", CONSTANT, MONOMIAL, &(this->get_region_ids()));
    system.add_variable("Vbot", CONSTANT, MONOMIAL, &(this->get_region_ids()));
  }
  else
  {
    system.add_variable("Vtop", FIRST, &(this->get_region_ids()));
    system.add_variable("Vbot", FIRST, &(this->get_region_ids()));
  }
  system.add_vector("currdens");
  system.attach_assemble_object(_my_assembly);
  system.init();

  // we use cm as length units here
  double mesh_units = 100 * get_mesh_units();
  get_scaling().set_calc_mesh_units(mesh_units);
}


void
PVModule::parse_options(void)
{

  string discr = get_option("discretization", "FEM");
  if (discr == "FEM")
    _discretization = FEM;
  else if (discr == "DEC")
    _discretization = DEC;
  else
    throw InitFailedException("Unknown discretization scheme: "
                              + discr);

  if (_discretization == DEC)
    throw InitFailedException("DEC discretization is not implemented yet.");


  _spice = get_option("spice_executable", _spice);
  get_parameter("voltage", _voltage);
 
}


void
PVModule::do_setup_solution_variables(void)
{
  // we declare our solution variables
  if (_discretization == DEC)
  {
    declare_solution(TopPotential, REAL, CELL, "V");
    declare_solution(BottomPotential, REAL, CELL, "V");
    declare_solution(CellPotential, REAL, CELL, "V");
    declare_solution(CurrentDensity, REAL, CELL, "A/cm^2");
    declare_solution(ContactCurrent, REAL, GLOBAL, "A");
  }
  else
  {
    declare_solution(TopPotential, REAL, NODES, "V");
    declare_solution(BottomPotential, REAL, NODES, "V");
    declare_solution(CellPotential, REAL, NODES, "V");
    declare_solution(CurrentDensity, REAL, NODES, "A/cm^2");
    declare_solution(ContactCurrent, REAL, GLOBAL, "A");
  }
}



void
PVModule::do_solve(void)
{

  TiberLinearSystem& sys = get_equation_system<TiberLinearSystem>();

  sys.assemble();

  // write Spice netlist
  string netlist = get_output_directory() + "/" + get_name() + "_spice.cir";
  ofstream of(netlist);
  of << "* PV Module scale simulation model \n" << std::endl;

  // idea: loop through matrix, and take coupling elements as resistors
  // DOF indices can be used as node names

  DofMap& dof_map =  sys.get_dof_map();
  const unsigned int vtop = sys.variable_number("Vtop");
  unsigned int n_dofs = dof_map.n_dofs();
  // we have two DoFs per node

  // the elementary cell subcircuit might add new nodes,
  // starting from this one
  unsigned int next_node_id = n_dofs + 1;

  vector<numeric_index_type> indices_t;
  vector<numeric_index_type> indices_b;
  vector<double> values_t;
  vector<double> values_b;

  // now we loop through the matrix rows. Conductance values that are
  // = 0 (no connection) can be ignored
  // First we connect all module nodes, then we add two additional nodes
  // for the voltage source and ground
  for (unsigned int i = 0; i < n_dofs; i += 2)
  {
    // top node id
    // increase by one because 0 is reserved for ground
    const unsigned int this_node = i + 1;

    // extract information from system matrix
    sys.matrix->get_row(i, indices_t, values_t);
    sys.matrix->get_row(i+1, indices_b, values_b);

    // indices_t.size() == indices_b.size() by construction
    for (unsigned int j = 0; j < indices_t.size(); ++j)
    {
      unsigned int other_node = indices_t[j] + 1;

      // diagonal element
      if (other_node == this_node)
      {
        // the area is given in cm^2
        double area = values_t[j];

        if (area > 0)
        {
          // create netlist of the elementary cell, scaled with area
          // bottom node
          unsigned int bot_node = this_node + 1;

          // for now we just take one of them, assuming that all
          // are equivalent
          PVModuleModel* pvm = dynamic_cast<PVModuleModel*>(
                                           *get_physical_models().begin());

          ElementaryCell* elementary = pvm->get_elementary_cell();

          // call the model
          // TODO: for now we do not have an element to pass
          // NOTE: next_node_id is expected to be a valid new node id
          // at return
          elementary->write_netlist(this_node, bot_node,
                                    next_node_id, area,
                                    nullptr,
                                    get_mesh().point(i/2), of);

        }

      }
      else if ((other_node > this_node) && (values_t[j] > 0.0))
      {
        // add top sheet resistance, or vertical R
        // this_node -- R -- other_node

        double resistance = 1.0 / values_t[j];   
        string t_or_v = "T ";

        if (other_node == (this_node + 1))
        {
          // P2 top-to-bottom connection
          t_or_v = "V ";
        }

        of << "R" << i / 2 << "_" << this_node << "_" << other_node
           << t_or_v << this_node << " " << other_node << " " << resistance << "\n";
      }


      // now the bottom layer
      unsigned int current_node = this_node + 1;
      other_node = indices_b[j] + 1;
      if ((other_node > this_node) && (values_b[j] > 0.0))
      {
        double resistance = 1.0 / values_b[j];   

        of << "R" << i / 2 << "_" << current_node << "_" << other_node
           << "B " << current_node << " " << other_node << " " << resistance << "\n";
      }

    }
  }

  // the node where voltage source is attached to
  unsigned int input_node = next_node_id;

  // now add a small resistor from voltage source node to each of _src_ids
  for (auto&& node : _src_ids)
  {
    of << "R" << input_node << "_" << node << " "
       << input_node << " " << node << " " << _rsource << "\n";
  }

  // and similarly from each of _gnd_id to "0"
  for (auto&& node : _gnd_ids)
  {
    of << "R" << "0" << "_" << node << " "
       << node << " 0 " << _rgnd << "\n";
  }
  
  of << "Vbias "<< input_node << " 0 DC " << _voltage << "\n";
  of << "* End of netlist \n";
  of << "* Simulation command \n";
  of << ".control \n";
  of << " op\n";
  //of << " dc Vbias 0 4 0.1 \n";

  string ngspice_res = get_output_directory() + "/" + get_name() + "_spice_output.txt";
  // this lets write the scale factor only once (e.g. the voltage in a spice sweep)
  of << "set wr_singlescale\n";
  of << " wrdata " << ngspice_res << " i(Vbias)";
  for (int nm = 1; nm <= n_dofs; nm++)
	    of << " V(" << nm << ")";
		  
		
  of << "\n";
  of << ".endc \n";
  of << ".end \n";
  of.flush();
  of.close();



  // call ngspice
  Messages::info("calling Spice: " + _spice);
  string log = get_output_directory() + "/" + get_name() + "_spice.log";
  string outfile = get_output_directory() + "/" + get_name() + "_spice.dat";
  //string cmdline = _spice + " -b -o " + log + " -r " + outfile + " " + netlist ;
  string cmdline = _spice + " -b -o " + log + " " + netlist + " > /dev/null";
  int ret = std::system(cmdline.c_str());

  if (ret == -1)
    throw(SolveFailedException("Could not run Spice."));

  // parse output and populate solution vectors
  Messages::info("parse ngspice results");

  // we need to put voltages into the solution vector
  libMesh::NumericVector<Number>& solution = sys.get_local_solution_vector();

  ifstream file(ngspice_res);
  string line;

  // in the current implementation we should have exactly one line
  while (std::getline(file, line))
  {
    stringstream ss(line);
    double value;
    int column_indx = 1;
    double src_vol;

    while (ss >> value)
    {
      if (column_indx == 2)        // column 2 represent the total current of the cell
        _current = -value; // we invert sign to get current referred to the simulated device 

      if (column_indx >= 3)
      {
        solution.set(column_indx - 3, value);
      }
      column_indx++;
    }

  }
  file.close();

  solution.close();
  sys.update();

  // now we can calculate the current density
  calculate_current_density();
  
  ostringstream os;
  os << "Voltage at terminal : " << _voltage << " V\n";
  os << "Current at terminal : " << _current << " A\n";
  os << "Module power        : " << _current*_voltage << " W\n";
  Messages m;
  m.indent();
  m.info(os.str());
  
}


void
PVModule::calculate_current_density(void)
{

  TiberLinearSystem& sys = get_equation_system<TiberLinearSystem>();
  const NumericVector<libMesh::Number>& solution = sys.get_solution_vector();
  NumericVector<libMesh::Number>& currdens = sys.get_vector("currdens");

  // loop through matrix, and take coupling elements as resistors
  // calculate Kirchhoff in each active node and devide by area

  DofMap& dof_map =  sys.get_dof_map();
  const unsigned int vtop = sys.variable_number("Vtop");
  unsigned int n_dofs = dof_map.n_dofs();

  vector<numeric_index_type> indices;
  vector<double> values;

  // we put the vertical current density on the first of the two DoFs
  for (unsigned int i = 0; i < n_dofs; i += 2)
  {
    // extract information from system matrix
    sys.matrix->get_row(i, indices, values);

    double voltage_i = solution(i);
    double current = 0.0;
    double area = 0.0;

    for (unsigned int jj = 0; jj < indices.size(); ++jj)
    {
      unsigned int j = indices[jj];

      if (j == i)
      {
        area  = values[jj];
      }
      else
      {
        double cond = values[jj];
        double voltage_j = solution(j);
        current += (voltage_i - voltage_j) * cond;
      }
    }

    // add current to voltage source
    if (_src_ids.count(i+1))
    {
      current += (voltage_i - _voltage) / _rsource;
    }

    // add current to ground 
    if (_gnd_ids.count(i+1))
    {
      current += voltage_i / _rgnd;
    }

    if (area > 0)
    {
      currdens.set(i, current / area);
    }
  }

  currdens.close();
}



void
PVModule::do_print_info(void)
{
  Messages::info("PVModule: lumped element model for photovoltaic module"
                 " simulation.");
}


PhysicalModel*
PVModule::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{
  return PVModuleModel::create(mat, options);
}



PhysicalModel*
PVModule::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  return PVModuleBoundaryModel::create(boundary, options);
}


void
PVModule::get_solution_secure(std::map<ID, std::vector<double> >& values)
{
  map<ID, vector<double> >::iterator mapit(values.begin());
  const map<ID, vector<double> >::iterator mapend(values.end());
  for ( ; mapit != mapend; ++mapit)
  {
    ID id = mapit->first;

    if (id == ContactCurrent)
    {
      values[id] = vector<double>(1, _current);
    }
  }

}


void
PVModule::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();
  const NumericVector<libMesh::Number>& currdens = system.get_vector("currdens");

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int v_top = system.variable_number("Vtop");
  const unsigned int v_bot = system.variable_number("Vbot");

  FEType fe_type = system.variable_type(v_top);
  unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dofs_top;
  vector<unsigned int> dofs_bot;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  //const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  //const vector<Point>& real_pts = fe->get_xyz();

  ID subdomain = elem->subdomain_id();

  fe->reinit(elem, &p);

  dof_map.dof_indices(elem, dofs_top, v_top);
  dof_map.dof_indices(elem, dofs_bot, v_bot);
  const unsigned int n_dofs = dofs_top.size();

  //PVModuleModel& mod = *get_bulk_model<PVModuleModel>(elem);


  for (unsigned int n = 0; n < np; n++)
  {
    double ut  = 0.0;
    double ub  = 0.0;
    double curr = 0.0;

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      ut += phi[i][n] * solution(dofs_top[i]);
      ub += phi[i][n] * solution(dofs_bot[i]);
      curr += phi[i][n] * currdens(dofs_top[i]);
    }

    if (values.count(TopPotential))
      values[TopPotential][n] = ut;

    if (values.count(BottomPotential))
      values[BottomPotential][n] = ub;

    if (values.count(CellPotential))
      values[CellPotential][n] = ut - ub;

    if (values.count(CurrentDensity))
      values[CurrentDensity][n] = curr;

  }
}


void
PVModule::plot_globaldata(void)
{
  /*
  string outdir = get_output_directory();
  if (!_currents.empty()){
    string filename(outdir + "/" + get_output_filename() + "_IV.dat");
    ofstream file;
    file.open(filename.c_str());
    if (file.good())
    {
      file << "# " << get_type() << " IV characteristic (" << get_name() << ")\n";
      file << "# " << "Voltage(V) Current(A)" << "\n";
      for (unsigned int i = 0; i < _voltages.size(); ++i)
      {
        file << _voltages[i] << " " << _currents[i] << endl; 
      }
      
	  file << "\n";

    }
    file.close();

  }
  */
	
}


void
PVModule::assemble(void)
{
  if (_discretization == DEC)
    assemble_dec_dual();
  else
    assemble_fem();
}



void
PVModule::assemble_dec_dual(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
 
  
  DofMap& dof_map =  system.get_dof_map();

  const unsigned int vtop = system.variable_number("Vtop");
  const unsigned int vbot = system.variable_number("Vbot");


}



void
PVModule::assemble_fem(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
 
  
  DofMap& dof_map =  system.get_dof_map();

  const unsigned int vtop = system.variable_number("Vtop");
  const unsigned int vbot = system.variable_number("Vbot");

  FEType fe_type = dof_map.variable_type(vtop);

  // the finite element
  unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));
  unique_ptr<QBase> qrule(QBase::build(QTRAP, dim));
  fe->attach_quadrature_rule(qrule.get());

  const vector<Real>& JxW = fe->get_JxW();
  const vector<Point>& q_point = fe->get_xyz();
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();

  /*
  // the surface finite element
  unique_ptr<FEBase> fe_face(build_finite_element(dim, fe_type));
  unique_ptr<QBase> qface(QBase::build(myopts.quadrature_type, dim-1, myopts.integration_order));
  fe_face->attach_quadrature_rule(qface.get());

  const vector<Real>& JxW_face = fe_face->get_JxW();
  const vector<Point>& qface_point = fe_face->get_xyz();
  const vector<vector<Real> >&  phi_face = fe_face->get_phi();
  const vector<vector<RealGradient> >& dphi_face = fe->get_dphi();
  const vector<Point>& normal = fe_face->get_normals();
  */

  vector<unsigned int> dof_indices;
  vector<unsigned int> dof_indices_top;
  vector<unsigned int> dof_indices_bot;

  DenseMatrix<Number> Ke;
  //DenseVector<Number> Fe;

  DenseSubMatrix<Number> Ket(Ke);


  MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

  for ( ; el != end_el ; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices(elem, dof_indices);
    const unsigned int n_dofs = dof_indices.size();

    dof_map.dof_indices(elem, dof_indices_top, vtop);
    dof_map.dof_indices(elem, dof_indices_bot, vbot);

    fe->reinit(elem);

    // resize the element matrix/rhs (does also zero them out)
    Ke.resize(n_dofs, n_dofs);

    Ket.reposition(0, 0, n_dofs/2, n_dofs/2);

    PVModuleModel& mod = *get_bulk_model<PVModuleModel>(elem);

    PVModuleModel::RegionType reg_type = mod.get_region_type();

    /* Idea:
     * the off-diagonal entries of the matrix will contain
     * the conductance values between different nodes (top, bottom,
     * and vertical connection). On the diagonal we put the area of
     * the active region of the node, used later for scaling the 
     * elementary cell parameters. 
     */

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {

      auto rsheet = mod.get_sheet_resistances(elem, q_point[qp]);
      double stop = 1.0 / rsheet.first;
      double sbot = 1.0 / rsheet.second;

      if (reg_type == PVModuleModel::P3)
        stop = 1e-9;

      if (reg_type == PVModuleModel::P1)
        sbot = 1e-9;

      double sconn = 1.0 / mod.get_connection_resistance(elem, q_point[qp]);

      for (unsigned int i = 0; i < n_dofs/2; i++)
      {
        // for the bottom layer
        unsigned int ii = i + n_dofs/2;

        for (unsigned int j = 0; j < n_dofs/2; j++)
        {
          if (j == i)
            continue;

          // top layer conductance contribution
          double topcon = JxW[qp] * (dphi[i][qp] * (stop * dphi[j][qp]));
          Ke(i, j) -= topcon;

          // bottom layer conductance contribution
          unsigned int jj = j + n_dofs/2;
          double botcon = JxW[qp] * (dphi[i][qp] * (sbot * dphi[j][qp]));
          Ke(ii, jj) -= botcon;
        }

        double area = JxW[qp] * phi[i][qp] * phi[i][qp];

        // active area contribution
        if (reg_type == PVModuleModel::ACTIVE)
          Ke(i, i) += area;

        // vertical connection conductance contribution
        if (reg_type == PVModuleModel::P2)
        {
          double vertconn = sconn * area;
          Ke(i, ii) += vertconn;
          Ke(ii, i) += vertconn;
        }
      }
    }


    // check the element sides, which may be associated with contacts
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      PVModuleBoundaryModel* mod_int =
        get_interface_model<PVModuleBoundaryModel>(elem, s);

      if (mod_int != NULL)
      {
        //auto side = elem->side_ptr(s);

        PVModuleBoundaryModel::ContactType type = mod_int->get_contact_type();
        PVModuleBoundaryModel::ContactLayer layer = mod_int->get_contact_layer();

        auto& id_set = (type == PVModuleBoundaryModel::GND) ? _gnd_ids : _src_ids;

        for (unsigned int n = 0; n < elem->n_nodes(); n++)
        {
          if (elem->is_node_on_side(n, s))
          {
            unsigned int dof_id = dof_indices_top[n];
            if (layer == PVModuleBoundaryModel::BOTTOM)
              dof_id = dof_indices_bot[n];
		  
            // spice node id 0 is reserved for ground
            dof_id++;

            id_set.insert(dof_id);

            if (layer == PVModuleBoundaryModel::BOTH)
              id_set.insert(dof_indices_bot[n]+1);
          }
        }
      }
    }
    

    //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
    system.matrix->add_matrix(Ke, dof_indices);
    //system.rhs->add_vector(Fe, dof_indices);

  }
  system.matrix->close();
  //system.matrix->print_matlab("K.m");

}
