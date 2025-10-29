// $Id: MyPoisson.C 2356 2011-02-19 22:54:42Z maufder $

#include "MolecularDynamics.h"
#include "TiberLinearSystem.h"
#include "Messages.h"

#include "MDModel.h"
#include "MDBoundaryModel.h"
#include "library.h"
#include "AtomisticStructure.h"
#include "SimulationEnvironment.h"

#include "equation_systems.h"
#include "dof_map.h"
#include "quadrature_gauss.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "tensor.h"
#include "fe_interface.h"
#include "Material.h"
#include "Database.h"
#include "Specie.h"
// This is needed in order to create the shared module library
// The first string is the class name of the object to be created,
// the second one is the name of the module as it should be referred
// in the input file (the Makefile defines MODULE_NAME, which can be used here).
TIBER_MODULE(MolecularDynamics, MODULE_NAME)


using namespace std;


MolecularDynamics*
MolecularDynamics::_this = NULL;


MolecularDynamics::MolecularDynamics(const ModelOptions& options) :
  SimulationInterface(options)
{
  // there's nothing to be done
}


MolecularDynamics::~MolecularDynamics(void)
{
  // there's nothing to be done
}


MolecularDynamics*
MolecularDynamics::create(const ModelOptions& options)
{
  // we could use the options to create different implementations
  // or something like that.
  return new MolecularDynamics(options);
}


void
MolecularDynamics::do_init(void)
{

  parse_options();
  create_equation_system("linear");

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  system.add_variable("T",FIRST);

  //system.attach_assemble_function(assemble);
  system.init();
  (system.solution)->add(0.0);

}


void
MolecularDynamics::parse_options(void)
{

  //Read lammps options
  myopts.potential_file = get_options().get_option("potential_file","pot_lammps");
  myopts.timestep = get_options().get_option("time_step",0.0);
  myopts.nvestep = get_options().get_option("nve_step",100);
  myopts.nvtstep = get_options().get_option("nvt_step",100);
  myopts.rescalestep = get_options().get_option("rescale_step",10);
  myopts.xbc = get_options().get_option("x_bc","p");
  myopts.ybc = get_options().get_option("y_bc","p");
  myopts.zbc = get_options().get_option("z_bc","p");
  myopts.Nevery = get_options().get_option("Nevery",10);
  myopts.Nrepeat = get_options().get_option("Nrepeat",10);
  myopts.Nfreq = get_options().get_option("Nfreq",100);

  myopts.dumpstep = get_options().get_option("dump_step",10);
  myopts.thermostep = get_options().get_option("thermo_step",10);

}


void
MolecularDynamics::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(LatticeTemp, REAL, NODES, "K");
  declare_solution(ThermalFlux, VECTOR, NODES, "W/cm^2");

}


void
MolecularDynamics::do_solve(void)
{
   double SimTemp = SimulationOptions::temperature;

  _this = this;
   vector<AtomisticStructure*> atom_structures;

   get_environment().get_device().get_atomistic_structures("all",atom_structures);

   atom_structures[0]->print_structure("before.xyz");

   double scale = atom_structures[0]->get_scale();

   void* lmp;

   lammps_open_no_mpi(0,NULL,&lmp);

   //CORE #1 -> Initialization
   cout<<lammps_command(lmp,(char*) "units metal")<<endl;
   cout<<lammps_command(lmp,(char*) "dimension 3")<<endl;

   //Boundary conditions
     {
       char buffer [100];
       int n=sprintf (buffer, "boundary %s %s %s",myopts.xbc.c_str(),myopts.ybc.c_str(),myopts.zbc.c_str());
       cout<<lammps_command(lmp,buffer)<<endl;
     }

   cout<<lammps_command(lmp,(char*) "boundary p p p")<<endl;
   cout<<lammps_command(lmp,(char*) "atom_style atomic")<<endl;
   cout<<lammps_command(lmp,(char*) "neighbor 2.0 bin")<<endl;
   cout<<lammps_command(lmp,(char*) "neigh_modify every 1")<<endl;

   //Time step
   {
     char buffer [100];
     int n=sprintf (buffer, "timestep %f",myopts.timestep);
     cout<<lammps_command(lmp,buffer)<<endl;
   }

   std::vector< Atom >& structure =  atom_structures[0]->get_structure_atoms();

   ID natoms = structure.size();

   //Initialize box
   double deltab = 1e-9;
   double box[6];
   box[0] = deltab;
   box[1] = -deltab;
   box[2] = deltab;
   box[3] = -deltab;
   box[4] = deltab;
   box[5] = -deltab;

   double positions[3*natoms];
   for (unsigned int na = 0; na < structure.size(); na++)
   {
     RealGradient pos = structure[na].get_position();

     if (pos(0)<box[0])
       box[0] = pos(0);

     if (pos(0)>box[1])
       box[1] = pos(0);

     if (pos(1)<box[2])
       box[2] = pos(1);

     if (pos(1)>box[3])
       box[3] = pos(1);

     if (pos(2)<box[4])
       box[4] = pos(2);

     if (pos(2)>box[5])
       box[5] = pos(2);

     positions[3*na]   = pos(0);
     positions[3*na+1] = pos(1);
     positions[3*na+2] = pos(2);

   }

   //Create box
   unsigned int n_types = atom_structures[0]->get_N_types();
   double delta = 1e-6;
   std::vector<double> masses(n_types);

   {
   char buffer [100];
   int n=sprintf(buffer, "region 1 block %f %f %f %f %f %f units box",box[0]-delta,box[1]+delta,box[2]-delta,box[3]+delta,box[4]-delta,box[5]+delta);
   cout<<lammps_command(lmp,buffer)<<endl;
   }
   {
   char buffer[100];
   int n=sprintf (buffer, "create_box %i 1",n_types);
   cout<<lammps_command(lmp,buffer)<<endl;
   }

   //CORE #2->Transfer atom coords to LAMMPS and write masses vector
   cout<<"Importing structure...";
   //cout<<"StructureSize: "<<atom_structures.size()<<endl;
   for (unsigned int na = 0; na < structure.size(); na++)
      {
        ID atom_type  = atom_structures[0]->get_type_index(structure[na].get_specie().get_string());
        masses[atom_type-1] = structure[na].get_specie().get_mass();
        char buffer [100];
        int n=sprintf (buffer, "create_atoms %i single %f %f %f units box",atom_type,positions[3*na],positions[3*na+1],positions[3*na+2]);
        char* str = lammps_command(lmp,buffer);
      }

   int natoms_lmp = lammps_get_natoms(lmp);
   cout<<"done. Number of atoms: "<<natoms_lmp<<endl;

   for (ID ns = 0;  ns < atom_structures[0]->get_N_types(); ns++)
   {
     char buffer [100];
     int n=sprintf (buffer, "mass %i %f",ns+1,masses[ns]);
     lammps_command(lmp,buffer);
   }

    //MODULE
    //get the potential file
    {
      std::string name_file = myopts.potential_file;
      char *a=new char[name_file.size()+1];
      a[name_file.size()]=0;
      memcpy(a,name_file.c_str(),name_file.size());
      lammps_file(lmp,a);
    }

    //Variables and computes
    lammps_command(lmp,(char*) "compute myKE all ke/atom");
    lammps_command(lmp,(char*) "variable mytemp atom c_myKE*2.0/3.0/8.6E-5");

    //This commands are needed in order to get the average temperature
    {
       char buffer [100];
       int n=sprintf (buffer,"fix AVETIME all ave/atom %i %i %i v_mytemp",myopts.Nevery,myopts.Nrepeat,myopts.Nfreq);
       lammps_command(lmp,buffer);

       lammps_command(lmp,(char*) "variable AVETEMP atom f_AVETIME");
    }

    //Compute variables
    if (myopts.nvtstep > 0)
    {
      {
         char buffer [100];
         int n=sprintf (buffer, "thermo %i",myopts.thermostep);
         lammps_command(lmp,buffer);
      }

      {
         char buffer [100];
         int n=sprintf (buffer, "velocity all create %T 429349 dist uniform ",SimTemp);
         //int n=sprintf (buffer, "velocity all create temp %T",SimTemp);
         lammps_command(lmp,buffer);
      }

      {
          char buffer [100];
          int n=sprintf (buffer, "fix NVT all nvt temp %f %f 0.0001 ",SimTemp,SimTemp);
          lammps_command(lmp,buffer);
      }

      {
        char buffer [100];
        int n=sprintf (buffer, "dump Dump all custom %i atom.lmp id type x y z v_mytemp v_AVETEMP",myopts.dumpstep);
        lammps_command(lmp,buffer);
      }

      {
        char buffer [100];
        int n=sprintf (buffer, "run %i",myopts.nvtstep);
        cout<<lammps_command(lmp,buffer)<<endl;
      }

    }

   //CORE #3 -> //Trasfer back coords to AtomisticStructure. NOTE: alla atom MUST be in the original box
   double coords[3*natoms_lmp];
   lammps_get_coords(lmp,coords);

   for (unsigned int na = 0; na < structure.size(); na++)
   {
     Tensor1 newpos(0);
     newpos(1) = coords[3*na];
     newpos(2) = coords[3*na+1];
     newpos(3) = coords[3*na+2];

     structure[na].set_position(newpos);
   }
   atom_structures[0]->print_structure("after.xyz");

   //CORE #4 -> get a general variable from LAMMPS and trasfer it to a nodal variable

    //if (myopts.nvestep + myopts.nvtstep > 0)
    if (myopts.nvtstep > 0)
    {


      double* lammps_temp = (double *) lammps_extract_variable(lmp, (char*) "AVETEMP", (char *) "all");
      //double * lammps_temp =(double *) lammps_extract_fix(lmp,(int) 5,(int) 5,(int) 5,(int) 5);
      //double* lammps_temp = (double *) lammps_extract_fix(lmp,5,1,0,0);
      int* id = (int *) lammps_extract_atom(lmp, (char*) "id");

      double temp[natoms];

      for (unsigned int na = 0; na < structure.size(); na++)
         {
           ID tag = id[na]-1;
           temp[tag] = lammps_temp[tag];
         }


    //write the temperature on nodes
    TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();
    system.set_options(get_solver_options());
    //const NumericVector<Number>& solution = system.get_solution_vector();

    NumericVector<Number>*  scaling = (system.solution)->clone().release();
    scaling->zero();

    const unsigned int dim = get_mesh().mesh_dimension();

    const DofMap& dof_map = system.get_dof_map();
    const unsigned int t_var = system.variable_number("T");

    FEType fe_type = system.variable_type(t_var);
    std::unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));

    vector<unsigned int> dof_indices;

    // element shape functions
    const vector<vector<Real> >& phi = fe->get_phi();
    const vector<Point>& real_pts = fe->get_xyz();
   //-----------------------------------------------------

    const unsigned int system_number = system.number();

    vector<Point> rel_point(1);
    for (unsigned int na = 0; na < structure.size(); na++)
    {

       const Elem* elem = structure[na].get_elem();
       vector<Point> tmp_point(1);
       tmp_point[0](0) = structure[na].get_position()(0)/scale;
       tmp_point[0](1) = structure[na].get_position()(1)/scale;
       tmp_point[0](2) = structure[na].get_position()(2)/scale;

       FEInterface::inverse_map(get_mesh().mesh_dimension(),FEType(),elem,tmp_point,rel_point);


       fe->reinit(elem, &rel_point);

       //cout<<endl;
       for (ID n = 0; n<elem->n_nodes(); n++ )
       {
         Node* node = elem->get_node(n);

         if (node->n_dofs(system_number, t_var) == 0)
            continue;

         ID n_dof = node->dof_number(system_number,t_var,0);

         double value = temp[na] * phi[n][0];

         //double value = phi[n][0]*300.0;
         system.solution->add(n_dof,value);

         //Build the scaling factor
         scaling->add(n_dof,phi[n][0]);
       }
      }

     //Scaling the shape function
     MeshBase::const_node_iterator  nd  = get_mesh().active_nodes_begin();
     const MeshBase::const_node_iterator nd_end = get_mesh().active_nodes_end();

     for ( ;  nd != nd_end ; ++nd)
     {
       Node* node = *nd;

       // If there are no DOFs, it's not a node of the simulation domain
       if (node->n_dofs(system_number, t_var) == 0)
         continue;

       ID n_dof = node->dof_number(system_number,t_var,0);

       double old_value = (*system.solution)(n_dof);

       double scale_factor = (*scaling)(n_dof);

       if (scale_factor > 0.0)
        system.solution->set(n_dof,old_value/scale_factor);
       else
        system.solution->set(n_dof,SimulationOptions::temperature);

     }
     //Release memory
     delete scaling;
    }

    //Delete the lammps object
    lammps_close(lmp);

}



void
MolecularDynamics::do_print_info(void)
{
  Messages::info("TiberCAD/LAMMPS interface");
}


PhysicalModel*
MolecularDynamics::create_bulk_model(const ModelOptions& options,
    const Material* mat) const
{

  return  MDModel::create(mat, options);

}

PhysicalModel*
MolecularDynamics::create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const
{
  return MDBoundaryModel::create(boundary, options);
}



void
MolecularDynamics::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{

  unsigned int np = p.size();

    TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

    const NumericVector<Number>& solution = system.get_solution_vector();

    const unsigned int dim = get_mesh().mesh_dimension();

    const DofMap& dof_map = system.get_dof_map();

    const unsigned int t_var = system.variable_number("T");

    FEType fe_type = system.variable_type(t_var);
    std::unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));

    vector<unsigned int> dof_indices;

    // element shape functions
    const vector<vector<Real> >& phi = fe->get_phi();
    const vector<vector<RealGradient> >& dphi = fe->get_dphi();
    const vector<Point>& real_pts = fe->get_xyz();

    ID subdomain = elem->subdomain_id();

    fe->reinit(elem, &p);

    dof_map.dof_indices(elem, dof_indices, t_var);
    const unsigned int n_dofs = dof_indices.size();

    MDModel& mod = *get_bulk_model<MDModel>(elem);

    bool do_temperature = values.count(LatticeTemp);


    for (unsigned int n = 0; n < np; n++)
    {
      Real T = 0.0;
      // do interpolation
      for (unsigned int i = 0; i < n_dofs; i++)
        T += phi[i][n] * solution(dof_indices[i]);

      if (do_temperature)
      {
          values[LatticeTemp][n] = T;
       }


    }


}



