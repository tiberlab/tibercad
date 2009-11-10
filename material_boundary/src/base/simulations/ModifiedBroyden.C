// $Id$

#include "ModifiedBroyden.h"
#include "newmatio.h"
#include "petsc.h"
using namespace std;
//---------------------------------------------------//

void ModifiedBroyden::do_solve(void)
{
 
  parse_options();

  clear_f();

  _it_number = 1;

  

  init_x();



  evaluate_and_save_f();
 
  

  calculate_kappa_matrix();
 
  calculate_new_first_iteraion_solution();


  for (   ; !_converged  ;   )
  {
    
    _it_number++;

    try
    {
     

      if (_it_number > get_maximum_iterations() ) 
        throw SolveFailedException("ModifiedBroyden::Number of iterations exceeded\n");
      
     
      do_step();
      
      
      double x = estimate_step();

      if (get_monitor())
      {

	double x = estimate_step();
	cout.flush();
	cout << "<<<<------------------------------------------------------------\n";
	cout << get_name() << " (Broyden-Johnson):  iteration " <<  _it_number 
	     <<"  Correction = " << x <<"  Error = " << _rel_error << "\n";
	cout << "--------------------------------------------------------------->>>>\n";
	cout.flush();
	
      }

      
	//double xd = _it_number;
	//double yd = log(_rel_error)/log(10);


    }      
    catch(RBD_COMMON::BaseException& e)
    {
      cerr << e.what() << "\n";
      throw SolveFailedException("Error in Mathematics");
    }

    
  }

  clear_f();

}

//------------------------------------------------------------//
inline void ModifiedBroyden::do_step(void)
{
  evaluate_and_save_f();
   
  calculate_kappa_matrix();

  update_omega_vector();

  calculate_mu_and_w_diag();
     
  calculate_tilde_a();
      
  calculate_lambda();

  update_tilde_a();
      
  calculate_a();
     
  calculate_beta_matrix();
       
  calculate_zeta_matrix();
   
  calculate_eta_vector();
     
  calculate_p_matrix();     
      
  calculate_new_solution();

}





//---------------------------------------------------------//

inline void ModifiedBroyden::calculate_tilde_a(void)
{
  NEWMAT::Matrix a_tilde_temp;

  if (_it_number > 2)
    a_tilde_temp = _a_tilde;

  _a_tilde.ReSize(_it_number - 1, _it_number - 1);

  if (_it_number > 2)
    _a_tilde.SymSubMatrix(1,_it_number-2) = a_tilde_temp;

  a_tilde_temp.Release();

  _a_tilde(_it_number - 1,_it_number - 1) = 1;

  for (unsigned int i = 1 ; i <= _it_number - 2; i++)
    _a_tilde(_it_number-1, i) = _lambda(i);
  
  
  

}
//------------------------------------------------------------//
inline void ModifiedBroyden::calculate_lambda(void)
{
  
  
  _lambda.ReSize(_it_number - 1);

  for (unsigned int i = 1; i <= _it_number - 1; i++)
    _lambda(i) = _mu(i,i) * (_kappa(_it_number, i+1) - _kappa(_it_number,i));
  

}
//------------------------------------------------------------//

inline void ModifiedBroyden::update_tilde_a(void)
{

  for (unsigned int i = 1; i <= _it_number - 2; i++)
    _a_tilde(_it_number - 1, i) = _mu(_it_number - 1, _it_number - 1)
      * (_lambda(i) -  _a_tilde(_it_number - 1, i));

  

}
//------------------------------------------------------------//
inline  void ModifiedBroyden:: calculate_a(void)
{
  _a.ReSize(_it_number - 1);

  for (unsigned int j = 1; j <= _it_number - 1; j++)
    for (unsigned int i = j; i <= _it_number - 1; i++)
      _a(i,j) = _w(i,i)*_w(j,j) * _a_tilde(i,j);
}
//------------------------------------------------------------// 

//------------------------------------------------------------//

inline void ModifiedBroyden::calculate_eta_vector(void)
{
  _eta = _w * _zeta * _lambda;
}

//------------------------------------------------------------//

inline void ModifiedBroyden::calculate_p_matrix(void)
{

  NEWMAT::Matrix _p_temp = _p;

  _p.ReSize(_it_number, _it_number);

  if (_it_number > 2)
    _p.SymSubMatrix(1, _it_number - 1) = _p_temp;

  _p_temp.Release();


  for (unsigned int i = 2 ; i <= _it_number; i++)
  {
    double t = 0;


    for (unsigned int j = i ; j <= _it_number - 1; j++)
    {
      t += _eta(j)  * _p (j,i);      
    }

   
    _p(_it_number,i) =  _eta(i - 1) - t;

  }
  

}

//------------------------------------------------------------------------------------//

inline void  ModifiedBroyden::calculate_new_solution(void)
{

  _x_correction->zero(); 

  _x_correction->add(_alpha * ( _p(_it_number, _it_number) - 1.0), *(_f[_it_number]) );

  for (unsigned int i = 2; i <= _it_number - 1; i++)
    _x_correction->add( _p(_it_number,i), *(_f[i]) );
     
  *_x += *_x_correction;

  
}

//-------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------//
inline double ModifiedBroyden::estimate_step()
{
  double t1 = 0;
  double t2 = 0;
  
  for (unsigned int i = 0; i < _vector_size; i++)
  {
    t1 += (*_x_correction)(i) * (*_x_correction)(i);

    t2 += (*_x)(i) * (*_x)(i) ;
  }
  
  t1 /= _vector_size;
  t2 /= _vector_size;

  

  return(sqrt(t1)/sqrt(t2));

}
//--------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------//
void ModifiedBroyden::parse_options(void)
{
 

  SelfconsistentSolver::parse_options();

  const ModelOptions& mod_opt = get_options();

  if (mod_opt.find_option("omega_0"))
  {
    _omega_0 = mod_opt.get_option("omega_0",1.0);
    _estimate_omega_0 = false;
  } 
  else
  {
    _estimate_omega_0 = true;
  }


  _alpha = mod_opt.get_option("relaxation_factor",0.3);


  _number_of_x_to_use = mod_opt.get_option("history_length",10);

  
  _scale_factor = mod_opt.get_option("scale_factor", 1);
  
  
}


//-------------------------------------------------------------------------------------//
void ModifiedBroyden::do_init(void)
{
  
  SelfconsistentSolver::do_init();


}
//-------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------//
inline void ModifiedBroyden::init_x(void)
{
  initialize();
  
  const NumericVector<double>& solution  = get_solution_vector(); 


  

  _vector_size = solution.size();

 
  if (_x.get() == NULL) 
  {
    _x = NumericVector<double>::build();
    _x->init(_vector_size);
  }
  else
  {
    _x->init(0);
    _x->init(_vector_size);
  }


  *_x = solution;

  _x->scale(_scale_factor);
 

  if (_x_correction.get() == NULL)
  {  
    _x_correction =  NumericVector<double>::build();
    _x_correction->init(_vector_size);
  }
  else
  {
    _x_correction->init(0);
    _x_correction->init(_vector_size);
  }
   

  
 
 


  if  (_estimate_omega_0)
    _omega_0 = _x->l2_norm();
 
}
