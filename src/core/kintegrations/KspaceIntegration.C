
// Basic include files needed for the mesh functionality.
#include "fe.h"
#include "fe_interface.h"
// Define generic quadrature rules.
#include "quadrature.h"
#include "quadrature_gauss.h"
#include "quadrature_trap.h"

// For mesh refinement
#include "mesh_refinement.h"
#include "kelly_error_estimator.h"
// Define the Finite Element object.

// Define useful datatypes for finite element
// matrix and vector components.
#include "KspaceIntegration.h"
#include "VTKIO.h"

#include "SimulationOptions.h"
#include "Messages.h"

using namespace std;

KspaceIntegration::KspaceIntegration(const ModelOptions& options)
 : TiberModelObject(options)
{
  _kspace = NULL;
  
}

//-------------------------------------------------------//
KspaceIntegration::~KspaceIntegration()
{
  delete _kspace;
}
//-------------------------------------------------------//
void KspaceIntegration::calculate_density()
{

  int verbose = get_option("verbose",0);

  Mesh* kmesh = const_cast <Mesh*>( _kspace->get_k_mesh() );
  //const Mesh* kmesh =  _kspace->get_k_mesh();


  unsigned int k_dim = kmesh->mesh_dimension();
    
  //-----------------------------------------------------------

  AutoPtr<FEBase> fe( FEBase::build(k_dim, FEType(fem_order) ));

  AutoPtr<QBase> qrule(QBase::build(quadrature_type, k_dim, integration_order));

  fe->attach_quadrature_rule(qrule.get());
cout<<"fe qrule attached\n";

  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<Point>& q_point = fe->get_xyz();


  MeshBase::const_element_iterator it_k_space= kmesh->active_elements_begin();
  const MeshBase::const_element_iterator it_k_end  = kmesh->active_elements_end();



  double factor = 1.0;
  {
    for (short i = 0; i < k_dim; i++)  factor /= (2.0 * M_PI);

    factor *= _kspace->get_degeneracy_factor() * opt.degeneracy;
  }

  //std::vector<unsigned int> dof_indices;
  // ------- INTEGRATION --------------------------------------------------------

  for ( ; it_k_space != it_k_end ; ++it_k_space) //loop over k space elements
  {
    const KElem* kelem = *it_k_space;

    KMeshToIntegratedValue::iterator it_k_elem;

    DofField dens_at_k_elem;
   
    // Skipping elements already computed in a previous mesh refin steps
    // kspace_integral == map<KMesh,double>;
    it_k_elem = kspace_integral.find(kelem);

    if (it_k_elem == kspace_integral.end())
    {
cout<<"fe reinit\n";
      fe->reinit(kelem);

cout<<"done\n";
      for (unsigned int qp=0; qp<q_point.size(); qp++)
      {//qp

	if (verbose > 2)
	{
          std::cout<<"(KIntegration) element "<< kelem->id()+1 <<"/" << kmesh->n_elem() <<
                        " point "<< qp+1 <<"/"<<q_point.size() << std::endl;
        }
	if (verbose > 3) 
		std::cout << "(KIntegration) k_point:  "<<  q_point[qp] << std::endl;

        DofField  dens_at_k_point;
        
        double integrated_quantity;

	calculate_for_k_point( q_point[qp], dens_at_k_point, integrated_quantity);


        // build the map between k-points and the integrated error quantity
	kspace_integral[kelem] += integrated_quantity * JxW[qp] * factor;

	// iterator over the real mesh 
	// add quad point contrib for every real-space element         
	DofField::iterator dens_it = dens_at_k_point.begin();
	DofField::iterator dens_it_end = dens_at_k_point.end();

	for ( ; dens_it != dens_it_end ; ++dens_it) //loop over real space elements
	{
	  const Elem* el = dens_it->first;

	  dens_at_k_elem[el] += (dens_it->second) * JxW[qp] * factor;
	} 

      }  //qp sum (dens_at_k_elem is computed)
 
      if (verbose > 3)
	  std::cout << "Contribution at k-element  "<<kspace_integral[kelem]<<"\n";

      // either register on a map (density_at_k) or update immediatly real_space_dens
      if (opt.k_domain_refinement)
	density_at_k.insert(pair<const KElem*, DofField >(kelem, dens_at_k_elem));
      else
      {
        DofField::iterator field_it  = dens_at_k_elem.begin();
        DofField::iterator field_end = dens_at_k_elem.end();
      
        for( ; field_it != field_end ; ++field_it )
          real_space_density[field_it->first]  += field_it->second;
       }

    } // if new k_elem
  
  } // end loop on active kelem

  //--------------------------------------------------------------------------//
  if (opt.k_domain_refinement)
  {
    MeshBase::const_element_iterator it_k_space = kmesh->active_elements_begin();
    const MeshBase::const_element_iterator it_k_end = kmesh->active_elements_end();

    for ( ; it_k_space != it_k_end ; ++it_k_space)
    {
      const KElem* kel = *it_k_space;

      DofField& dens_at_k_elem = density_at_k[kel];

      DofField::iterator field_it =   dens_at_k_elem.begin();
      DofField::iterator field_end =  dens_at_k_elem.end();

      for ( ; field_it != field_end ; ++field_it )
 	real_space_density[field_it->first] += field_it->second;
    }

  }

}


//-------------------------------------------------------------------------------//
void KspaceIntegration::calculate_convergent_density()
{

  int verbose = SimulationOptions::verbose();

  Mesh* kmesh = const_cast <Mesh*> (_kspace->get_k_mesh() );

  cout <<"(KIntegration) Calculate k-integral "<<endl;

  density_at_k.clear();

  kspace_integral.clear();

  real_space_density.clear();

  calculate_density();


  if (opt.k_domain_refinement)
  {

    if (verbose > 1)
      std::cout << "Simulation " << get_name() << " " << "is performing k space refinement\n";

    MeshRefinement mesh_refinement(*kmesh);

    double norm_of_error = opt.relative_accuracy + 1.0;

    while(norm_of_error > opt.relative_accuracy)
    {
      if (opt.uniform_refinement)
	mesh_refinement.uniformly_refine(1);
      else
      {

	ErrorVector error = ErrorVector(kmesh->n_elem(), kmesh);

	estimate_error_for_refinement(error);

	mesh_refinement.refine_fraction() = opt.refine_fraction;

	mesh_refinement.max_h_level() = opt.maximum_ref_level;

	mesh_refinement.coarsen_fraction() = 0.0;

	if (verbose > 3)
	{
	  cout << "\nError vector:  " << error.size() << "\n";
	  for (unsigned int i = 0; i < error.size(); i++ )
	    cout << setprecision(10) <<  error[i] << "\n";
	}

	mesh_refinement.flag_elements_by_error_fraction(error);

	mesh_refinement.refine_and_coarsen_elements();

	if (verbose > 1)
	{
	  cout << "\nThe new mesh: \n";
	  kmesh->print_info();
	}

	old_density = real_space_density;

	real_space_density.clear();

	calculate_density();

	norm_of_error = estimate_error();

	if (verbose > 1)
	  std::cout <<  "\n\n Relative Error: " << norm_of_error<<"\n"<<std::endl;

      } // h-refinement

    } //refinement loop

  }//end of refinement block

}

//------------------------------------------------------------------------------------//

void KspaceIntegration::parse_options( )
{

  const ModelOptions& mod_opt = get_options();

  fem_order = FIRST;

  string quad_type = mod_opt.get_option("quadrature_type","gaussian");

  if(quad_type == "gaussian") quadrature_type = QGAUSS;
  else if(quad_type == "trapezoidal") quadrature_type = QTRAP;
  else throw  InitFailedException("Kspace: unsupported quadrature type: "+quad_type+"\n" ); 

  string int_order = mod_opt.get_option("quadrature_order","third");

  if(int_order == "first") integration_order = FIRST;
  else if(int_order == "second") integration_order = SECOND;
  else if(int_order == "third") integration_order = THIRD;
  else if(int_order == "fourth") integration_order = FOURTH;
  else if(int_order == "fifth") integration_order = FIFTH;
  else if(int_order == "sixth") integration_order = SIXTH;
  else if(int_order == "seventh") integration_order = SEVENTH;
  else if(int_order == "eighth") integration_order = EIGHTH;
  else throw  InitFailedException("Kspace: unsupported quadrature order: "+int_order+"\n" );   
   
  opt.uniform_refinement      = mod_opt.get_option("uniform_refinement",false);

  opt.refine_fraction         = mod_opt.get_option("refine_fraction", 0.3);
  opt.maximum_ref_level       = mod_opt.get_option("maximum_ref_level", 8);
  opt.relative_accuracy       = mod_opt.get_option("relative_accuracy", 1e-2);

  opt.degeneracy                = mod_opt.get_option("degeneracy",1);
  opt.k_domain_refinement       = mod_opt.get_option("refine_k_space", false);
  opt.log_output                = mod_opt.get_option("log_output",  false);

  //additional_name_suffix  = mod_opt.get_option("suffix", "");

  //cout<<"done"<<endl;

}
//------------------------------------------------------------------------------------//

void KspaceIntegration::do_solve( )
{

  calculate_convergent_density();
}

//--------------------------------------------------------------------------------------//
void KspaceIntegration::do_init(void)
{
  
  std::cout<<"(KSI) k-int: "<<std::endl;
  ModelOptions kopts;

  if(has_option("mesh_units"))
   kopts.set_option("mesh_units",get_option("mesh_units",0.0));
  else  
   throw InitFailedException("K-integration internal error: mesh_units must be initialized");	  

  if(has_option("k_space_dimension"))
   kopts.set_option("k_space_dimension",get_option("k_space_dimension",0));
  else
   throw InitFailedException("K-integration internal error: k_space_dimension must be initialized");	  

  std::vector<unsigned int>  num_nodes;
 
  if (has_option("number_of_nodes"))
  {
     get_option("number_of_nodes",num_nodes);
     kopts.set_option("number_of_nodes", num_nodes);
  }
  else if (has_option("number_of_elements"))
  {
     get_option("number_of_elements",num_nodes);
     for(int i=0; i< num_nodes.size(); i++)
           if(num_nodes[i]>0) ++num_nodes[i];

     kopts.set_option("number_of_nodes", num_nodes);
  }
  else
   throw InitFailedException("K-integration internal error: number_of_nodes must be initialized");	  

  if (has_option("wedge"))
    kopts.set_option("wedge", get_option("wedge",""));

  
  kopts.set_option("k_space_basis", get_option("k_space_basis",true));


  double k_max = get_option("k_max",0.1);
  std::cout<<"(KSI) k_max: "<< k_max<<std::endl;

  kopts.set_option("k_max",k_max);


  std::vector<double> k_vector(3,0.0);   
  k_vector[0]=0.0;     k_vector[1]=0.0;     k_vector[2]=k_max; 

  get_option("k1", k_vector);
  kopts.set_option("k1",k_vector);

  k_vector[0]=0.0;     k_vector[1]=k_max;     k_vector[2]=0.0; 
  get_option("k2", k_vector);
  kopts.set_option("k2",k_vector);

  k_vector[0]=k_max;     k_vector[1]=0.0;     k_vector[2]=0.0; 
  get_option("k3", k_vector);
  kopts.set_option("k3",k_vector);  
    
  kopts.set_option("mesh_order",get_option("mesh_order","first"));

  std::cout<<"(KSI) kspace init: "<<std::endl;

  _kspace = new Kspace(kopts);

  if(_kspace==NULL)
    throw InitFailedException("Could not initialize k-space");
  else
    Messages::info("k-space initialized");


  parse_options();

}

//--------------------------------------------------------------------------------------//
void KspaceIntegration::estimate_error_for_refinement(ErrorVector& error) 
{

  std::fill(error.begin(), error.end(), 0.0);

  const Mesh* kmesh = _kspace->get_k_mesh();

  MeshBase::const_element_iterator       elem_it1  = kmesh->elements_begin();
  const MeshBase::const_element_iterator elem_end1 = kmesh->elements_end();

  for (; elem_it1 != elem_end1; ++elem_it1)
  {
    const KElem* el = *elem_it1;
    const unsigned int el_id = el->id();

    error[el_id] = abs(kspace_integral[el]); //test
  }

}

//---------------------------------------------------------------------------------
double  KspaceIntegration::estimate_error(void) 
{

  double result;
  double t1 = 0.0; double t2 = 0.0;

  DofField::iterator it = real_space_density.begin();
  const DofField::iterator it2 = real_space_density.end();

  for (; it != it2; ++it)
  {
    const Elem* el = it->first;

    t1 += real_space_density[el] * real_space_density[el];

    t2 += (real_space_density[el] - old_density[el]) *  
          (real_space_density[el] - old_density[el]);

  }

  result = t2/t1;

  return(result);

}
//--------------------------------------------------------------------------------------//
unsigned int KspaceIntegration::count_elements() const
{
  unsigned int result = 0;

  const Mesh* kmesh = _kspace->get_k_mesh();


  MeshBase::const_element_iterator       elem_it  = kmesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = kmesh->active_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    result++;

  return result;

}
//=================================================================//
/*std::vector<double>   KspaceIntegration::get_density_in_k_space(void)  const
/{


  const Mesh* kmesh = _kspace->get_k_mesh();


  MeshBase::const_element_iterator       elem_it  = kmesh->active_elements_begin();
  const MeshBase::const_element_iterator elem_end = kmesh->active_elements_end();

  const unsigned int n_active_elements = count_elements();

  vector<double> result( n_active_elements );


  unsigned int j = 0;
  for ( ; elem_it !=  elem_end; ++elem_it )
  {
    const KElem* el = *elem_it;


    KMeshToIntegratedValue::const_iterator it1 = kspace_integral.find(el);


    result[j] = it1->second / el->volume();


    j++;
  }


  return(result);

}
*/





