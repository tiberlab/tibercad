#include "ModifiedBroyden.h"
using namespace std;
//---------------------------------------------------//

void ModifiedBroyden::do_solve(void)
{
/*
  _it_number = 1;

  init_X();

  do_iteration();

  get_X_from_iteration();

  unsigned int n = _X.size(); 

 
  evaluate_F();
 

  calculate_kappa_matrix();

  calculate_new_solution();

  pass_X_to_iteration();


  bool converged = false;


  for (   ; !converged  ;   )
  {
    do_iteration();

    get_X_from_iteration();

    _vector_size = _X.size();

    evaluate_F();

    calculate_kappa_matrix();
    
    calculate_mu_and_w_diag();

    calculate_tilde_a();

    calculate_lambda();

    calculate_a();

    calcvulate_beta_matrix();

    calculate_zeta_matrix();

    calculate_eta_vector();

    calculate_p_matrix();

    calculate_new_solution();

    pass_X_to_iteration();
    

    double x = evaluate_difference();

    if (x < get_relative_tolerance() ) converged = true;

  }
*/
}


//------------------------------------------------------------//

inline void ModifiedBroyden::evaluate_F(void)
{
  vector<double> temp(_vector_size);
  for (unsigned int i = 0; i < _vector_size; i++ )
    temp[i] = _X[i] - _X1[i];

  _F.push_back(temp);
  
  
}

//------------------------------------------------------------//

inline void ModifiedBroyden::calculate_kappa_matrix(void)
{
/*
  TNT::Fortran_Array2D<double> k_temp;
  if (_it_number > 1)  k_temp = _kappa.copy();

  _kappa(_it_number, _it_number);
  
  
  for (unsigned int i = 1; i <= _it_number; i++)
  {
    double t = 0;
    for (unsigned int j = 0; j <  _vector_size; j++)
      t += _F[i][j] * _F[_it_number][j];
     _kappa(_it_number,  i) = t;
     _kappa( i, _it_number) = t;

  }

  if (_it_number > 1) 
    for (unsigned int i = 1; i <= _it_number; i++)
      for (unsigned int j = 1; j <= _it_number; j++)
	_kappa(i,j) =  k_temp(i,j);

*/
}
//--------------------------------------------------------//
inline void ModifiedBroyden::calculate_mu_and_w_diag(void)
{
/*
  //--------------------------------------------------------
  //calculation of mu
  TNT::Fortran_Array1D<double> mu_temp;

  mu_temp = _mu.copy();


  _mu = TNT::Fortran_Array1D<double>(_it_number - 1);

  _mu(_it_number - 1) = 1.0/sqrt(
				_kappa(_it_number,_it_number) - 
				2.0*_kappa(_it_number,_it_number-1) +
				_kappa(_it_number-1,_it_number-1) 
				);

  for (unsigned int i = 1; i <= _it_number - 2; i++)
    _mu(i) = mu_temp(i);



  //calculation of w
 
  _w = TNT::Fortran_Array1D<double>(_it_number - 1, 1.0);
*/  
}
//--------------------------------------------------------//
