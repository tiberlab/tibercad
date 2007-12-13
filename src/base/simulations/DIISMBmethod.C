#include "DIISMBmethod.h"

#include "newmatio.h"
inline
void DIISMBmethod::calculate_new_solution()
{
  _X_correction->zero(); 

  _X_correction->add(_alpha * ( _s(_it_number, _it_number) - 1.0), *(_F[_it_number]) );

  for (unsigned int i = 2; i <= _it_number - 1; i++)
  {
    _X_correction->add( _s(_it_number,i), *(_F[i]) );
  
    cerr << _s(_it_number,i) << "\n";
  }


  // _X_correction->add(- _alpha  , *(_F[_it_number]) );
   
  *_X += *_X_correction;
}


inline 
void DIISMBmethod::calculate_tilde_a(void)
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
    _a_tilde(_it_number-1, i) = _mu(_it_number - 1, _it_number - 1)*( _lambda(_it_number,i) - _lambda(_it_number - 1,i)) ;

}


inline 
void DIISMBmethod::calculate_lambda_matrix(void)
{

  NEWMAT::Matrix lambda_temp;

  if (_it_number > 2)
    lambda_temp = _lambda;

  
  _lambda.ReSize(_it_number, _it_number - 1);

  
   if (_it_number > 2)
     _lambda.SubMatrix(1,_it_number-1, 1,  _it_number - 2 ) = lambda_temp;

   
   lambda_temp.Release();

   for (unsigned int i = 1 ; i <= _it_number - 1; i++)
   {
     _lambda(i, _it_number - 1) = _mu( _it_number - 1, _it_number - 1 ) * ( _kappa( _it_number,i) - _kappa( _it_number - 1,i) );
     _lambda(_it_number, i) = _mu(i,i) * ( _kappa( _it_number,i + 1) - _kappa( _it_number,i));
   }
   

}


inline
void DIISMBmethod::calculate_eta_matrix(void)
{
  _eta = _lambda * _w * _zeta;
}


inline 
void DIISMBmethod::calculate_q_matrix(void)
{
  NEWMAT::Matrix q_temp;
  
  if (_it_number > 2)
    q_temp = _q;

  _q.ReSize(_it_number , _it_number );


   if (_it_number > 2)
     _q.SymSubMatrix(1,_it_number - 1) = q_temp;

   q_temp.Release();

   for (unsigned int i = 1 ; i <= _it_number - 1; i++)
   {
     double t = 0;
     for (unsigned int j = 2 ; j <= _it_number - 1; j++)
       t += _s( _it_number - 1,j) * _kappa(j, _it_number);

     _q(_it_number - 1,i) = t;


     t = 0;

     for (unsigned int j = 2 ; j <= i; j++)
       t += _s( i,j) * _kappa(j, _it_number);

     _q(i, _it_number) = t;
   }
}

inline 
void DIISMBmethod::calculate_sigma_matrix(void)
{
 NEWMAT::Matrix sigma_temp;

 if (_it_number > 2)
    sigma_temp = _sigma;

 _sigma.ReSize(_it_number - 1 , _it_number);

 if (_it_number > 2)
   _sigma.SubMatrix(1,_it_number - 2, 1, _it_number - 1) = sigma_temp;
 
 sigma_temp.Release();

 for (unsigned int i = 1 ; i <= _it_number - 1; i++)
 {
   _sigma(_it_number - 1, i) = _alpha * (_q(_it_number - 1,i) - _kappa(_it_number,i));
   _sigma(i,_it_number) = _alpha * (_q(i,_it_number) - _kappa(i+1,_it_number));
 }
}

inline
void DIISMBmethod::calculate_nu_matrix(void)
{
  NEWMAT::SymmetricMatrix nu_temp;
  if (_it_number > 2)
    nu_temp = _nu;

  _nu.ReSize(_it_number - 1);

  if (_it_number > 2)
    _nu.SymSubMatrix(1,_it_number-2) = nu_temp;

  nu_temp.Release();

  for (unsigned int i = 1 ; i <= _it_number - 1; i++)
  {
    double t = 0;
    for (unsigned int j = 2 ; j <= i; j++)
      t += _s(i,j)*_q(_it_number - 1,j);

    _nu(i, _it_number - 1) = _alpha * _alpha * (_kappa(_it_number, i + 1) - _q(_it_number  - 1, i+1) - _q(i,_it_number) + t);
  }

  
}

inline 
void DIISMBmethod::calculate_DIIS_matrix(void)
{
 
  
  if (_it_number > _number_of_x_to_use)
    _eps.ReSize( _number_of_x_to_use, _number_of_x_to_use);
  else
    _eps.ReSize(_it_number, _it_number );

  

  if (_it_number > _number_of_x_to_use )
    _eps = (_alpha*_alpha) * _kappa.SymSubMatrix(_it_number - _number_of_x_to_use + 1,  _it_number);
  else
  {
    
    _eps = (_alpha*_alpha) * _kappa;

  }

  if (_it_number > _number_of_x_to_use )
    _eta_m = _eta.Rows(_it_number - _number_of_x_to_use + 1, _it_number);
  else
    _eta_m = _eta;


  

  NEWMAT::Matrix _sigma_m;

  if (_it_number > _number_of_x_to_use )
    _sigma_m =  _sigma.Columns(_it_number - _number_of_x_to_use + 1, _it_number);
  else
    _sigma_m = _sigma;

 

  _eps += _alpha * ( _eta_m * _sigma_m  + (_eta_m * _sigma_m).t());

 

  _eps += _eta_m * _nu * _eta_m.t();

 

  _sigma_m.Release();

}

inline 
void  DIISMBmethod::calculate_c_m(void)
{
  unsigned int size;

  if (_number_of_x_to_use > _it_number)
    size = _it_number;
  else
    size = _number_of_x_to_use;

  NEWMAT::Matrix _A(size + 1,size + 1);
  
  _A.SymSubMatrix(2, size + 1) = _eps;

  _A.SubMatrix(1,1,  2, size + 1) = -1;
 
  _A.SubMatrix(2, size + 1, 1, 1) = -1;


  cerr << setprecision(18) << _eps << "\n";

  NEWMAT::ColumnVector _b(size + 1);
 
  _b(1) = -1;
  
  NEWMAT::ColumnVector _qq(size + 1);
  _qq = _A.i()*_b; 

  cerr << "_qq\n";
  cerr << _qq << "\n";

  _c_m = _qq.Rows(2,size + 1);

  _A.Release();
  _b.Release();
  _qq.Release();
}


inline
void  DIISMBmethod::calculate_o_m_and_phi(void)
{

  unsigned int size;

  if (_number_of_x_to_use > _it_number)
    size = _it_number;
  else
    size = _number_of_x_to_use;

  _o_m.ReSize(_it_number - 1, size );



  _o_m = 0;

  for (unsigned int i = 1 ;  i <= size ;  i++)
  {

    cerr << _it_number - 1 << "   " <<  _it_number - 1 - (size - i) + 1 <<"\n";
    for (unsigned int j = _it_number - 1 ; 
	 (j >= _it_number - 1 - (size - i) + 1) && (j >= 1)  ;
	 j--)
    {
     
     
      _o_m(j,i) = 1;
    
    }
  }


  

  _phi = (_eta_m.t() + _o_m)*_c_m;

 
}


inline 
void DIISMBmethod::calculate_s_matrix(void)
{

  NEWMAT::Matrix _s_temp = _s;

  _s.ReSize(_it_number, _it_number);

  if (_it_number > 2)
    _s.SymSubMatrix(1, _it_number - 1) = _s_temp;

  _s_temp.Release();


  for (unsigned int i = 2 ; i <= _it_number; i++)
  {
    double t = 0;


    for (unsigned int j = i ; j <= _it_number - 1; j++)
    {
      t += _phi(j)  * _s(j,i);      
    }

   
    _s(_it_number,i) =  _phi(i - 1) - t;

  }
  

 


}

inline
void DIISMBmethod::do_step(void)
{
 
  cerr <<"1\n";
  evaluate_and_save_F();
 

  calculate_kappa_matrix();
 

  update_omega_vector();
 

  calculate_mu_and_w_diag();
  

  calculate_lambda_matrix();

  cerr <<"2\n";
  calculate_tilde_a();


  calculate_a();
 

  calculate_beta_matrix();
  
  calculate_zeta_matrix();

  calculate_eta_matrix();
 
  cerr <<"3\n";
  calculate_q_matrix();
 

  calculate_sigma_matrix();
 

  calculate_nu_matrix();
  cerr <<"3.5\n";

  calculate_DIIS_matrix();
 
  cerr <<"3.8\n";
  calculate_c_m();

  cerr << "4\n";
  calculate_o_m_and_phi();
  cerr << "4.1\n";
  calculate_s_matrix();
  cerr << "4.2\n";
  calculate_new_solution();

}
