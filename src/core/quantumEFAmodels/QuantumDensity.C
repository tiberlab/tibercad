
#include "QuantumDensity.h"


using namespace std;

//============================================//

QuantumDensity::QuantumDensity()
{

}

//============================================//

void QuantumDensity::set_options( QuantumDensity::options& options   )
{
  opt =  options;
}


//============================================//

QuantumDensity:: QuantumDensity( EnvelopFunctionApprox* model )
{
  quantum_model = model;
  real_space_density_size = quantum_model->get_number_of_active_cells();
}

//============================================//

QuantumDensity:: QuantumDensity( EnvelopFunctionApprox* model, QuantumDensity::options& options  )
{
  quantum_model = model;
  opt = options;
  real_space_density_size = quantum_model->get_number_of_active_cells();

}



//============================================//
QuantumDensity:: ~QuantumDensity()
{
  delete kmesh;
}


//============================================//
std::vector<double>& QuantumDensity::get_density(void)
{
  return(real_space_density);
}

//============================================//
void QuantumDensity::build_k_grid()
{


  //build mesh
  kmesh = new Mesh(k_dim);

  ElemType type;

  if (k_dim == 1)
    {
      type = EDGE3;
    }

  if (k_dim == 2)
    {
      type = QUAD8;
    }

  if (k_dim == 3)
    {
      type = HEX27;
    }


  MeshTools::Generation::build_cube (*kmesh, 
				     num_nodes[0], num_nodes[1], num_nodes[2], 
				     kmin[0], kmax[0], 
				     kmin[1], kmax[1], 
				     kmin[2], kmax[2],
				     type);
  cerr << kmin[0] <<"   "<<  kmax[0] <<"   "<<  kmin[1] <<"   "<<  kmax[1] << "  "<<  kmin[2] <<"   "<<
    kmax[2] << "  " <<  num_nodes[0] <<"   "<<  num_nodes[1] << "   "  << num_nodes[2] << "\n";


  
  rotate_mesh(kmesh, transform_matrix);

  kmesh->print_info();


  

}
//


//============================================//
void QuantumDensity::calculate_density()
{

 


  //----------------------------------------------------------------------//
  //build equation system object

  EquationSystems eq(*kmesh);
      
  eq.add_system<ImplicitSystem> ("k-integration");
  
  

  ImplicitSystem& system = eq.get_system<ImplicitSystem> ("k-integration");

  system.add_variable("u", SECOND);

  const DofMap& dof_map = system.get_dof_map();
    
  FEType fe_type = dof_map.variable_type(0);
    
  AutoPtr<FEBase> fe (FEBase::build(k_dim, fe_type));
    
  QGauss qrule (k_dim, FIRST);
    
  fe->attach_quadrature_rule (&qrule);

  
  const std::vector<Real>& JxW = fe->get_JxW();
  
  const std::vector<Point>& q_point = fe->get_xyz();
  
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  MeshBase::element_iterator       it     = kmesh->active_elements_begin();
  const MeshBase::element_iterator it_el  = kmesh->active_elements_end();


  real_space_density.resize(real_space_density_size,0.0);

  for ( ; it != it_el ; ++it) 
    {

      const Elem* elem = *it;
      fe->reinit (elem);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
	{//qp
	  unsigned int num_nodes =  phi.size();
	  for (unsigned int i = 0; i <  num_nodes; i++)
	    {
	      const Node* nd = elem->get_node(i);
	      map <const Node*, vector <double> > :: iterator node_it;

	      node_it = k_point_density.find(nd);
	      if (node_it != k_point_density.end())
		{
		  vector<double>&  dens_at_k_node = node_it->second;

		  for (unsigned int j = 0; j <  real_space_density_size; j++)
		    {
		      real_space_density[j] += phi[i][qp] * dens_at_k_node[j];
		    }
		    

		}
	      else
		{
		  cerr << "WARNING! QuantumDensity   node is missing\n";
		}
	      

	    }

	}
      
    }


  //--------------------------------------------------------------------------//

  double pi = 4.0 * atan(1);

  double factor = 1;

  for (short i = 0; i < k_dim; i++)  factor /= (2.0 * pi);

  for (unsigned int j = 0; j <  real_space_density_size; j++)  real_space_density[j] *= factor;
   


  //--------------------------------------------------------------------------//

}
//================================================================//
void  QuantumDensity::define_k_space(Tensor1 k_vector, unsigned int n)
{

  double norm_k = norm(k_vector);
  kmin[0] = -norm_k/2.0;  kmax[0] = norm_k/2.0; num_nodes[0] = n;

  

  for (short i = 1; i < 4; i++)
    {
      transform_matrix(i,1) = k_vector(i)/norm_k; 
    } 

  k_dim = 1;
}

//================================================================//
void QuantumDensity::define_k_space(Tensor1 k_vector1,unsigned int n,  Tensor1 k_vector2, unsigned int m)
{

  double norm_k1 = norm(k_vector1);
  kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;


  
  double norm_k2 = norm(k_vector2);
  kmin[1] = -norm_k2/2.0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;

  
  for (short i = 1; i < 4; i++)
    {
      transform_matrix(i,1) = k_vector1(i)/norm_k1; 
      transform_matrix(i,2) = k_vector2(i)/norm_k2;
    }

  k_dim = 2;

}

//================================================================//

void QuantumDensity::define_k_space(Tensor1 k_vector1, unsigned int n, Tensor1 k_vector2, 
				    unsigned int m, Tensor1 k_vector3, unsigned int k)
{

  double norm_k1 = norm(k_vector1);
  kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;


  
  double norm_k2 = norm(k_vector2);
  kmin[1] = -norm_k2/2.0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;


  
  double norm_k3 = norm(k_vector3);
  kmin[2] = -norm_k3/2.0;  kmax[2] = norm_k3/2.0; num_nodes[2] = k;


  for (short i = 1; i < 4; i++)
    {
      transform_matrix(i,1) = k_vector1(i)/norm_k1; 
      transform_matrix(i,2) = k_vector2(i)/norm_k2;
      transform_matrix(i,3) = k_vector3(i)/norm_k3;
    }

  k_dim = 3;

}

//================================================================//
void QuantumDensity::calculate_at_each_k_point()
{

  build_k_grid();

  MeshBase::node_iterator       it     = kmesh->active_nodes_begin();
  const MeshBase::node_iterator it_el  = kmesh->active_nodes_end();

  double band_edge = quantum_model->get_band_edge();




  

 

  for ( ; it != it_el ; ++it) 
    {
      
      const Node*  nd  = *it;

      map<const Node*, vector<double> >::iterator it1;

      it1 =  k_point_density.find(nd);

      if (it1 ==  k_point_density.end() )

	{
	  double k_vector[3];


	  k_vector[0] = (*nd)(0);
	  k_vector[1] = (*nd)(1);
	  k_vector[2] = (*nd)(2);
	  
	  std::cout <<  k_vector[0] << "   " <<  k_vector[1] << "    " <<  k_vector[2] << "\n";
	  
	  quantum_model->apply_k_vector(k_vector);
	  
	  vector<double>  dens = quantum_model->calculate_convergent_density(opt.Temperature);


     
     
	  k_point_density.insert( pair<const Node*, vector<double> > (nd, dens) );
	}

    }


  

}


//==============================================================//

void QuantumDensity::rotate_mesh(Mesh* mesh, Tensor2Gen& RotMatrix)
{

  Tensor1 vec1;

  for (unsigned int n=0; n < mesh->n_nodes(); n++)
    {
       const Point p = mesh->node(n);
  
       vec1(1) = p(0);
       vec1(2) = p(1);
       vec1(3) = p(2);

       vec1 = RotMatrix * vec1;

       mesh->node(n) = Point( vec1(1), vec1(2), vec1(3) );

     }

}




//=============================================================//
