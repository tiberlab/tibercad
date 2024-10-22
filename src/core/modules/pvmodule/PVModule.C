// $Id: PVModule.C 5514 2024-05-03 15:22:47Z maufder $

#include "PVModule.h"
#include "PVModuleModel.h"
#include "PVModuleBoundaryModel.h"
#include "TiberLinearSystem.h"
#include "Messages.h"
#include "Database.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dense_matrix.h"
//#include "libmesh/dense_vector.h"
//#include "libmesh/dense_submatrix.h"
//#include "libmesh/dense_subvector.h"
//#include "libmesh/vector_value.h"

// This is needed in order to create the shared module library
#include "TiberModule.h"
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
  system.add_variable("Vtop", FIRST);
  system.add_variable("Vbot", FIRST);
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

  _spice = get_option("spice_executable", _spice);
  //_spice = get_option("Voltage", _voltage);
 
 
  // reading jv_ref.dat file, as reference for JV curve
  Database db;
  db.set_material("jv_ref", get_option("jv_ref", ""));
  ifstream ifs;
  ifs.open(db.get_data_file().c_str());
  if (ifs.fail() || !ifs.good())
    throw InitFailedException("Cannot read spectrum "
        "from file " + db.get_data_file());

  size_t i = 0;
  const size_t buf_len = 256;
  char buf[buf_len];
   while (ifs.good())
   {
     if (i == _jv_ref_v.size())
     {
       size_t n_new = _jv_ref_v.size() + 10;
       _jv_ref_v.reserve(n_new);
       _jv_ref_j.reserve(n_new);
     }
     ifs.getline(buf, buf_len);
     if (buf[0] != '#')
     {
       istringstream in(buf);
       double l, s;
       if (in >> l >> s)
       {
         _jv_ref_v.push_back(s);
         _jv_ref_j.push_back(l);
         i++;
       }
     }
   }
   ifs.close();
}


void
PVModule::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(TopPotential, REAL, NODES, "V");
  declare_solution(BottomPotential, REAL, NODES, "V");
  declare_solution(CurrentDensity, REAL, NODES, "A/cm^2");
}

void
PVModule::check_contact_node(unsigned int &id) const
{
  if (_gnd_ids.count(id))
    id = *_gnd_ids.begin();
  else if (_src_ids.count(id))
    id = *_src_ids.begin();
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
  // DOF indices can be used as node numbers

  DofMap& dof_map =  sys.get_dof_map();
  const unsigned int vtop = sys.variable_number("Vtop");
  unsigned int n_dofs = dof_map.n_dofs();
  // we have two DoFs per node

  // the elementary cell subcircuit will add new nodes,
  // starting from this one
  unsigned int next_node_id = n_dofs;

  vector<numeric_index_type> indices;
  vector<double> values;

  // new we loop through the matrix rows. Conductance values that are
  // < 1e-7 (no connection) can be ignored
  for (unsigned int i = 0; i < n_dofs; i += 2)
  {
    // extract information from system matrix
    sys.matrix->get_row(i, indices, values);

    for (unsigned int j = 0; j < indices.size(); ++j)
    {
      unsigned int other_node = indices[j];

      if (other_node == i)
      {
        // the area is given in cm^2
        double area = values[j];

        if (area > 0)
        {
          // create netlist of the elementary cell, scaled with area
          // increase by one because 0 is reserved for ground
          unsigned int this_node = i + 1;
          unsigned int bot_node = i + 2;

          // substitute node ids, if they are on contacts
          check_contact_node(this_node);
          check_contact_node(bot_node);

          // Defining a voltage-dependent current source based on the
          // JV-Ref file and the area of the element

          if (this_node != bot_node)
          {
            of << "B" << i / 2 << " " << this_node << " " << bot_node
               << " I=pwl(V(" << this_node << ")-V(" << bot_node << ")";

            for (int nm = 0; nm < _jv_ref_v.size(); nm++)
              of << ", " << _jv_ref_j[nm] << ", " << _jv_ref_v[nm] * area << "m";

            of << ")\n"
               << std::endl;
          }

          // adjust next_node_id
        }

      }
      else if (values[j] > 1e-7)
      {
        // add bottom or top sheet resistance
        // this_node -- R -- other_node

        double resistance = 1.0 / values[j];   
		    unsigned int this_node = i;
        string t_or_b = "T ";

        if (other_node < i)
        {
          // bottom
          this_node++;
          other_node++;
          t_or_b = "B ";
        }

        // increase by one because 0 is reserved for ground
        this_node++;
        other_node++;

        // substitute node ids, if they are on contacts
        check_contact_node(this_node);
        check_contact_node(other_node);

        if (this_node != other_node)
          of << "R" << i / 2 << "_" << this_node << "_" << other_node
             << t_or_b << this_node << " " << other_node << " " << resistance << "\n";
      }

    }
  }
  
  //of << "Vbias "<<*_gnd_ids.begin() <<" " << *_src_ids.begin() << " " << _voltage <<" \n";
  of << "Vbias "<<*_gnd_ids.begin() <<" " << *_src_ids.begin() << " " << "0" <<" \n";
  of << "* End of netlist \n";
  of << "* Simulation command \n";
  of << ".control \n";
  of << "	op\n";
  of << "	dc Vbias 0 4 0.1 \n";

  string ngspice_res = get_output_directory() + "/" + get_name() + "_spice_output.txt";
  of << "	wrdata " << ngspice_res << " i(Vbias)";
  for (int nm = 1; nm <= n_dofs; nm++ )
	  if (!_gnd_ids.count(nm) && !_src_ids.count(nm))
	    of << " V(" << nm << ")";
		  
		
  of << "\n";
  of << ".endc \n";
  of << ".end \n";
  of.close();




  for (auto&& a : _gnd_ids)
    std::cerr  << a << " ";
  std::cerr << "\n\n";

  for (auto&& a : _src_ids)
    std::cerr  << a << " ";
  std::cerr << "\n";

  // call ngspice
  Messages::info("calling Spice: " + _spice);
  string log = get_output_directory() + "/" + get_name() + "_spice.log";
  string outfile = get_output_directory() + "/" + get_name() + "_spice.dat";
  //string cmdline = _spice + " -b -o " + log + " -r " + outfile + " " + netlist ;
  string cmdline = _spice + " -b -o " + log + " " + netlist ;
  int ret = std::system(cmdline.c_str());

  if (ret == -1)
    throw(SolveFailedException("Could not run Spice."));

  // parse output and populate solution vectors
  Messages::info("parse output");
  ifstream file(ngspice_res);
  cout<<ngspice_res<<endl;
  string line;
  if (std::getline(file, line)) { // for test,  just reading the first line, voltage == 0 
    stringstream ss(line);
    double value;
	int column_indx=1;
	double src_vol;
	
	// parse output for given voltage and populate in an queue
	queue<double> temp_res;
	while (ss >> value) {
	  if (column_indx == 4)  // column 4 represent the total current of the cell
		  _current.push_back(value * 1e7); //mA/cm^2
	  if (column_indx == 5) // column 5 represent source voltage
	    src_vol = value;
	  if (column_indx > 5 && column_indx % 2 == 0){  // ignoring odd column, they are just repeating the source voltage 
	    temp_res.push(value);
	  }
      column_indx++;
    
	}
	// map the ngspice results to actual elements
	for (int nm = 1; nm <= n_dofs; nm++)
	  if (_gnd_ids.count(nm))
	    _spic_res.push_back(0);
	  else if (_src_ids.count(nm))
        _spic_res.push_back(src_vol);
	  else{
	    _spic_res.push_back(temp_res.front()); //results
	    temp_res.pop();
	}		
			
  }
  file.close();
  
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
PVModule::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();

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

    // do interpolation
    for (unsigned int i = 0; i < n_dofs; i++)
    {
      ut += phi[i][n] * solution(dofs_top[i]);
      ub += phi[i][n] * solution(dofs_bot[i]);
    }

    if (values.count(TopPotential))
      values[TopPotential][n] = ut;

    if (values.count(BottomPotential))
      values[BottomPotential][n] = ub;

  }
}


void
PVModule::plot_globaldata(void)
{
  string outdir = get_output_directory();
  if (!_current.empty()){
    string filename(outdir + "/" + get_output_filename() + "_CurrentDensity.dat");
    ofstream file;
    file.open(filename.c_str());
    if (file.good())
    {
      file << "# " << get_type() << " CurrentDensity (" << get_name() << ")\n";
      file << "# " << 1 << " CurrentDensity" << "\n";
      file << "# " << "CurrentDensity" << "\n";
      for (auto i :_current)
        file << i << endl; 
      
	  file << "\n";

    }
    file.close();

  }
	
}


void
PVModule::assemble(void)
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


  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

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
    //Fe.resize(n_dofs);

    Ket.reposition(0, 0, n_dofs/2, n_dofs/2);

    PVModuleModel& mod = *get_bulk_model<PVModuleModel>(elem);

    PVModuleModel::RegionType reg_type = mod.get_region_type();

    /* Idea:
     * the off-diagonal entries of the matrix will contain
     * the conductance values between different nodes (top, bottom,
     * and vertical connection). On the diagonal we put the area of
     * the active region of the node, used later for scaling the 
     * elementary cell parameters. Top and vertical connections are
     * put with j > i, bottom ones with j < i
     */

    // loop over the quadrature points
    for (unsigned int qp = 0; qp < qrule->n_points(); qp++)
    {

      auto rsheet = mod.get_sheet_resistances(elem, q_point[qp]);
      double stop = 1.0 / rsheet.first;
      double sbot = 1.0 / rsheet.second;

      if (reg_type == PVModuleModel::P3)
        stop = 0.0;

      if (reg_type == PVModuleModel::P1)
        sbot = 0.0;

      double sconn = 1.0 / mod.get_connection_resistance(elem, q_point[qp]);

      for (unsigned int i = 0; i < n_dofs/2; i++)
      {
        for (unsigned int j = i + 1; j < n_dofs/2; j++)
        {
          unsigned int ii = i;
          unsigned int jj = j;
          if (dof_indices[j] < dof_indices[i])
          {
            ii = j;
            jj = i;
          }

          // top layer conductance contribution
          Ke(ii, jj) -= JxW[qp] * (dphi[i][qp] * (stop * dphi[j][qp]));

          // bottom layer conductance contribution
          // Ke(i+n_dofs/2, j+n_dofs/2) += JxW[qp] * (dphi[i][qp] * (sbot * dphi[j][qp]));
          Ke(jj, ii) -= JxW[qp] * (dphi[i][qp] * (sbot * dphi[j][qp]));
          }

          // active area contribution
          if (reg_type == PVModuleModel::ACTIVE)
            Ke(i, i) += JxW[qp] * phi[i][qp] * phi[i][qp];

          // vertical connection conductance contribution
          if (reg_type == PVModuleModel::P2)
            Ke(i, i + n_dofs / 2) += JxW[qp] * sconn * phi[i][qp] * phi[i][qp];
      }
    }


    // the sides
    for (unsigned int s = 0; s < elem->n_sides(); s++)
    {
      PVModuleBoundaryModel* mod_int =
        get_interface_model<PVModuleBoundaryModel>(elem, s);

      if (mod_int != NULL)
      {
        //auto side = elem->side_ptr(s);

        PVModuleBoundaryModel::ContactType type = mod_int->get_contact_type();
        PVModuleBoundaryModel::ContactLayer layer = mod_int->get_contact_layer();

        for (unsigned int n = 0; n < elem->n_nodes(); n++)
        {
          if (elem->is_node_on_side(n, s))
          {
            unsigned int dof_id = dof_indices_top[n];
            if (layer == PVModuleBoundaryModel::BOTTOM)
              dof_id = dof_indices_bot[n];

            // spice node id 0 is reserved for ground
            dof_id++;
            
            if (type == PVModuleBoundaryModel::GND)
              _gnd_ids.insert(dof_id);
            else
              _src_ids.insert(dof_id);
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
