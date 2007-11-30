// $Id$

#include "ModifiedBroyden.h"
using namespace std;
//---------------------------------------------------//

void ModifiedBroyden::do_solve(void)
{

  _it_number = 1;

  init_X();

  do_iteration();

  get_X_from_iteration();


  evaluate_and_save_F();
 

  calculate_kappa_matrix();

  calculate_new_first_iteraion_solution();

  


  bool converged = false;


  for (   ; !converged  ;   )
  {
    do_iteration();

    get_X_from_iteration();

    _vector_size = _X.size();

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
    

    double x = estimate_error();

    if (x < get_relative_tolerance() ) converged = true;

  }

}


//------------------------------------------------------------//

inline void ModifiedBroyden::evaluate_and_save_F(void)
{
  vector<double> temp(_vector_size);
  for (unsigned int i = 0; i < _vector_size; i++ )
    temp[i] = _X[i] - _X1[i];

  _F.insert( pair<unsigned int, vector<double> > (_it_number, temp) );
  
  
}

//------------------------------------------------------------//

inline void ModifiedBroyden::calculate_kappa_matrix(void)
{

  NEWMAT::SymmetricMatrix k_temp;

  if (_it_number > 1)  k_temp = _kappa;

  _kappa.ReSize(_it_number);

  _kappa.Inject(k_temp);

   k_temp.Release();
    
   for (unsigned int i = 1; i <= _it_number; i++)
   {
     double t = 0;
     
     for (unsigned int j = 0; j <  _vector_size; j++)
       t += _F[i][j] * _F[_it_number][j];
     
     _kappa(_it_number,  i) = t;
     

   }


}
//--------------------------------------------------------//
inline void ModifiedBroyden::calculate_mu_and_w_diag(void)
{

  //--------------------------------------------------------
  //calculation of mu
  NEWMAT::DiagonalMatrix mu_temp;

  mu_temp = _mu;

  _mu.ReSize(_it_number - 1);

  _mu.Inject(mu_temp);

  mu_temp.Release();

  _mu(_it_number - 1) = 1.0/sqrt(
				 _kappa(_it_number,_it_number) - 
				 2.0*_kappa(_it_number,_it_number-1) +
				 _kappa(_it_number-1,_it_number-1) 
				 );

 

  //--------------------------------------------------------//
  //calculation of w 
  NEWMAT::DiagonalMatrix w_temp;

  w_temp = _w;

  _w.ReSize(_it_number - 1);

  _w.Inject(w_temp);

  w_temp.Release();
  
  _w(_it_number - 1) =  _omega(_it_number - 1);

  
  
}
//---------------------------------------------------------//

inline void ModifiedBroyden::calculate_tilde_a(void)
{
  NEWMAT::Matrix a_tilde_temp;

  a_tilde_temp = _a_tilde;

  _a_tilde.ReSize(_it_number - 1, _it_number - 1);

  _a_tilde.Inject(a_tilde_temp);

  a_tilde_temp.Release();

  _a_tilde(_it_number - 1,_it_number - 1) = 1;

  for (unsigned int i = 1 ; i <= _it_number - 2; i++)
    _a_tilde(_it_number-1, i) = _lambda(i);
  
  
  

}
//------------------------------------------------------------//
inline void ModifiedBroyden::calculate_lambda(void)
{
  
  
  _lambda.ReSize(_it_number - 1);

  for (unsigned i = 1; i <= _it_number - 1; i++)
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
   NEWMAT::SymmetricMatrix t(_it_number -1);
  
  t = (_omega_0 * I + _a);

  _beta = t.i();

}
//------------------------------------------------------------//
inline void ModifiedBroyden::calculate_zeta_matrix(void)
{

  NEWMAT::Matrix _beta_primed = _beta.SubMatrix(1,_it_number-1,1,_it_number-2);

  _zeta = _beta * _w * _mu + (_omega_0 * _omega_0) * _beta_primed * _zeta;
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

  _p.Inject(_p_temp);

  _p_temp.Release();

  for (unsigned int i = 1 ; i < _it_number; i++)
  {
    double t = 0;

    for (unsigned int j = 1 ; j < _it_number - 1; j++)
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
      t += _p(_it_number,i)*_F[i][k];
    

    _X_correction[k] = _alpha * (_p(_it_number, _it_number) - 1.0)*_F[_it_number][k] + _alpha*t;

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
inline double ModifiedBroyden::estimate_error()
{
  double t1 = 0;
  double t2 = 0;
  
  for (unsigned int i = 0; i < _vector_size; i++)
  {
    t1 += _X_correction[i] * _X_correction[i];
    t2 += _X[i] * _X[i];
  }
  
  return(t1/t2);

}
//--------------------------------------------------------------------------------------//
inline void ModifiedBroyden::update_omega_vector(void)
{
  _omega.ReSize(_it_number - 1);

  _omega = 1;

  for (unsigned int i = 1 ; i <= _it_number -  _number_of_x_to_use ; i++)
  {
    _omega(i) = 0;
  }

}
//--------------------------------------------------------------------------------------//
