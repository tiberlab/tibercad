// $Id$

#include "ModifiedBroyden.h"
#include "newmatio.h"
using namespace std;
//---------------------------------------------------//

void ModifiedBroyden::do_solve(void)
{
 

  parse_options();


  _it_number = 1;

  

  init_X();

 

  

  evaluate_and_save_F();
 
  

  calculate_kappa_matrix();
 
  calculate_new_first_iteraion_solution();

  
 



  for (   ; !_converged  ;   )
  {
    
    _it_number++;

    try{
     

      if (_it_number > get_maximum_iterations() ) 
	throw SolveFailedException("ModifiedBroyden::Number of iterations exceeded\n");
      
     
    

      evaluate_and_save_F();

   
      calculate_kappa_matrix();

    
    
     

      update_omega_vector();

     
      calculate_mu_and_w_diag();
     
      calculate_tilde_a();
      
      calculate_lambda();

      update_tilde_a();
      
      
    

      calculate_a();
     
    
      calcvulate_beta_matrix();
      
    
      calculate_zeta_matrix();

    
     
      calculate_eta_vector();

   
     
      calculate_p_matrix();

     
      
      calculate_new_solution();
      
      
      double x = estimate_step();

      if (get_monitor())
      {
	double x = estimate_step();
	cout.flush();
	cout << "<<<<------------------------------------------------------------\n";
	cout << get_name() << " (Broiden-Johnson):  iteration " <<  _it_number 
	     <<"  Correction = " << x <<"  Error = " << _rel_error << "\n";
	cout << "--------------------------------------------------------------->>>>\n";
	cout.flush();
      }

      

    }      
    catch(RBD_COMMON::BaseException& e)
    {
      cerr << e.what() << "\n";
      throw SolveFailedException("Error in Mathematics");
    }
    

  }

}


//------------------------------------------------------------//

inline void ModifiedBroyden::evaluate_and_save_F(void)
{
  


  NumericVector< double > & solution_before	= get_solution_vector ();


  //-------------//
  //set solution from previuos Broyden step

  for (unsigned int i = 0; i < _vector_size; i++)
    solution_before.set(i, _X[i]);

  //--------------//
  //do calculation

  solve_simulations();


  NumericVector< double > & solution_after	= get_solution_vector ();



  vector<double> temp(_vector_size);

  for (unsigned int i = 0; i < _vector_size; i++ )
  {
    temp[i] = _X[i] - solution_after(i);  
  }

  //----------------------
  //check if converged
  double t1 = 0;
  double t2 = 0;
  
  for (unsigned int i = 0; i < _vector_size; i++)
  {
    t1 += temp[i] * temp[i];

    t2 += _X[i] * _X[i] ;
  }


  _rel_error = sqrt(t1/t2);

  if (_rel_error <=  get_relative_tolerance())
    _converged = true;
  else
    _converged = false;

 
  //----------------------


  _F.insert( pair<unsigned int, vector<double> > (_it_number, temp) );
  
  
}

//------------------------------------------------------------//

inline void ModifiedBroyden::calculate_kappa_matrix(void)
{



  NEWMAT::SymmetricMatrix k_temp;

  if (_it_number > 1)  k_temp = _kappa;

  _kappa.ReSize(_it_number);

  if (_it_number > 1) 
  {
   
    _kappa.SymSubMatrix(1,_it_number - 1) = k_temp;

   
  }
     
  k_temp.Release();


  for (unsigned int i = 1; i <= _it_number; i++)
  {
    double t = 0;
     
    for (unsigned int j = 0; j <  _vector_size; j++)
    {
      
      t +=  _F[i][j] * _F[_it_number][j];
    }
    _kappa(_it_number,  i) = t;
     

  }
  
 
 
}
//--------------------------------------------------------//
inline void ModifiedBroyden::calculate_mu_and_w_diag(void)
{

  //--------------------------------------------------------
  //calculation of mu

  
  _mu.ReSize(_it_number - 1);

  for (unsigned int i = 1; i <= _it_number - 1; i++)
  {
     _mu(i, i) = 1.0/sqrt(	 _kappa(i+1,i+1) -  2.0*_kappa(i+1,i) +  _kappa(i,i) );

     _mu(i,i) = sqrt(_mu(i,i));//I think it should be, but not present in the paper!
  }


  //--------------------------------------------------------//
  //calculation of w 
  NEWMAT::DiagonalMatrix w_temp;

  if (_it_number > 2)
    w_temp = _w;

  _w.ReSize(_it_number - 1);

  if (_it_number > 2) 
    _w.SymSubMatrix(1,_it_number-2) = w_temp;

  w_temp.Release();
  
  _w(_it_number - 1) =  _omega(_it_number - 1);


  
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
inline void ModifiedBroyden::calcvulate_beta_matrix(void)
{
   NEWMAT::DiagonalMatrix I(_it_number -1);

   I = 1;

   NEWMAT::SymmetricMatrix t(_it_number -1);
  
   t = (_omega_0 * I + _a);


  

  _beta = t.i();

}
//------------------------------------------------------------//
inline void ModifiedBroyden::calculate_zeta_matrix(void)
{

  NEWMAT::Matrix _zeta_old = _zeta;

  _zeta.ReSize(_it_number - 1, _it_number - 1);

  _zeta = _beta * _w * _mu;

  if (_it_number > 2) 
  {
    NEWMAT::Matrix _beta_primed = _beta.SubMatrix(1,_it_number-1,1,_it_number-2);
    
    _zeta.SubMatrix(1,_it_number-1,1,_it_number-2) += (_omega_0 * _omega_0) * _beta_primed * _zeta_old;

    _beta_primed.Release();
  }

  _zeta_old.Release();


}
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

  for (unsigned int k = 0; k < _vector_size; k++)
  {
    double t = 0;

    for (unsigned int i = 2; i <= _it_number - 1; i++)
    {
      t += _p(_it_number,i)*_F[i][k];
     
    }

    //-----------------------------------------------
    //Broyden-Johnson step

     _X_correction[k] = _alpha * ( _p(_it_number, _it_number) - 1.0) * _F[_it_number][k] + _alpha*t; 

     //--------------------------------------------

     //-------------------------------------------
     //test -- will be linear relaxation 
     //_X_correction[k] =  _alpha * (- 1.0)*_F[_it_number][k]  ;  
     //------------------------------------------
    _X[k] +=  _X_correction[k];

   

  }

  

}

//-------------------------------------------------------------------------------------//
inline void ModifiedBroyden::calculate_new_first_iteraion_solution(void)
{
 for (unsigned int k = 0; k < _vector_size; k++)
 {
   _X[k] -= _alpha * _F[1][k];
 }
  
}
//-------------------------------------------------------------------------------------//
inline double ModifiedBroyden::estimate_step()
{
  double t1 = 0;
  double t2 = 0;
  
  for (unsigned int i = 0; i < _vector_size; i++)
  {
    t1 += _X_correction[i] * _X_correction[i];

    t2 += _X[i] * _X[i] ;
  }
  
  t1 /= _vector_size;
  t2 /= _vector_size;

  

  return(sqrt(t1)/sqrt(t2));

}
//--------------------------------------------------------------------------------------//
inline void ModifiedBroyden::update_omega_vector(void)
{
 

  _omega.ReSize(_it_number - 1);

  _omega = 1;


  if (_it_number > _number_of_x_to_use)
  {
    int number_of_zeros =  _it_number -  _number_of_x_to_use; 

    for ( int i = 1 ;  i <= number_of_zeros ; i++)
      _omega(i) = 0;
    
  }

  
}
//--------------------------------------------------------------------------------------//
void ModifiedBroyden::parse_options(void)
{
 

  SelfconsistentSolver::parse_options();

  const ModelOptions& mod_opt = get_options();

  _omega_0 = mod_opt.get_option("omega_0",1.0);

  _alpha = mod_opt.get_option("alpha",0.15);

  _number_of_x_to_use = mod_opt.get_option("history_length",10);

  
  
  
  
}


//-------------------------------------------------------------------------------------//
void ModifiedBroyden::do_init(void)
{
  
  SelfconsistentSolver::do_init();


}
//-------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------//
inline void ModifiedBroyden::init_X(void)
{
  initialize();
  
  const NumericVector<double>& solution  = get_solution_vector(); //after iteration //should be




  _vector_size = solution.size();

  _X.resize(_vector_size);

  for (unsigned int i = 0; i < _vector_size; i++ )
  {
      _X[i] = solution(i);
   

  }

  _X_correction.resize(_vector_size);
  
}
