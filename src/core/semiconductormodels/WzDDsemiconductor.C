using namespace std;
#include "WzDDsemiconductor.h"
const double WzDDsemiconductor::Hartree;
//--------------------------------------------------//
WzDDsemiconductor::WzDDsemiconductor(void): DDsemiconductor()
{
  par.Ev = 0;
  par.EgGamma = 0;
  par.m_c_zz = 0; //conduction band mass along [0001] direction
  par.m_c_xx = 0; //conduction band mass along [1-210] direction
  par.A1 = 0;
  par.A2 = 0;
  par.A4 = 0;
  par.A5 =0 ;
  par.A6=0;
  par.A7=0; //k.p parameters
  par.a_x=0;
  par.a_z =0; //hydrostatic deformation potential
  par.D1 = 0;
  par.D2 = 0;
  par.D3 = 0;
  par.D4 = 0;
  par.D5 =0 ;
  par.D6 =0; //deformation potentials
  par.delta_s =0;
  par.delta_cr =0; //sp

}

//--------------------------------------------------//

WzDDsemiconductor::WzDDsemiconductor(const WzDDparameters& params): DDsemiconductor()
{
  par = params;
}
//---------------------------------------------------//
void WzDDsemiconductor::set_Eg(const double EgGamma1)
{
  par.EgGamma = EgGamma1;
}
//---------------------------------------------------//
void WzDDsemiconductor::set_mass_Gamma(const double m_c_zz_, const double m_c_xx_)
{
  par.m_c_zz = m_c_zz_;
  par.m_c_xx = m_c_xx_;
}
//----------------------------------------------------//
void WzDDsemiconductor::set_deform_pot_cond(const double a_x_, const double a_z_)
{
  par.a_x = a_x_;
  par.a_z = a_z_;
}
//----------------------------------------------------/
void WzDDsemiconductor::calculate_conduction_band_extremum(void)
{
  vector<DDsemiconductor::band_extremum> result;
  
  double Ev_top;
  //---------------------------------
  //valence top reference energy
  //-------------------------------
  double d1 =  par.delta_cr;
  double d2 =  par.delta_s;
  double d3 = d2;
  
  double E1 = d1 + d2;

  double E2 = (d1 - d2)/2.0 + sqrt( (d1-  d2/2.0)*( d1- d2/2.0) + 2.0 * d3 * d3 );

  if (E1 > E2)
    Ev_top = par.Ev + E1;
  else
    Ev_top = par.Ev + E2;

  //------------------------------

  double energy = Ev_top + par.EgGamma;
 
  if (strained) energy += (strain(1,1) + strain(2,2))* par.a_x + par.a_z *  strain(3,3);
   
  DDsemiconductor::band_extremum band_ext;
  band_ext.energy = energy;
  band_ext.degeneracy = 2;
  band_ext.mass_DOS = pow(par.m_c_xx * par.m_c_xx * par.m_c_zz, 1/3) ;

  result.push_back(band_ext);

  conduction_band = result;
}
//----------------------------------------------------/

void  WzDDsemiconductor::calculate_valence_band_extremum(void)
{
  vector<DDsemiconductor::band_extremum>   result;

  DDsemiconductor::band_extremum extremum;

  vector<Tensor1> k_vector;
  //-----------------------------------------------
  Tensor1 k;
  k(1) = 0 ; k(2) = 0; k(3) = 0;
  k_vector.push_back(k);

  k(1) = k_max ; k(2) = 0; k(3) = 0;
  k_vector.push_back(k);
 
  k(1) = 0 ; k(2) = k_max; k(3) = 0;
  k_vector.push_back(k);

  k(1) = 0 ; k(2) = 0; k(3) = k_max;
  k_vector.push_back(k);

  k(1) =k_max  ; k(2) = k_max; k(3) = 0;
  k_vector.push_back(k);


  k(1) =k_max  ; k(2) = 0; k(3) = k_max;
  k_vector.push_back(k);

  k(1) =0.0    ; k(2) = k_max; k(3) = k_max;
  k_vector.push_back(k);


  //--------------------------------------------------
  vector< vector<double> >  eigenvalue = calculate_vb_bulk_states(k_vector);
    

  Tensor2Sym imass;
  for (short ind = 0; ind < 6; ind++)
    {
      if (eigenvalue[0][ind] + energy_cutoff > eigenvalue[0][5]) 
	{
	  extremum.degeneracy = 1;

	  extremum.energy  =eigenvalue[0][ind] ;
	  
	  imass(1,1) = (2.0 *(eigenvalue[0][ind] - eigenvalue[1][ind] )) / Hartree /(k_max * k_max);
	 
	 
	  imass(2,2) = (2.0 *(eigenvalue[0][ind] - eigenvalue[2][ind] )) / Hartree /(k_max * k_max);
	 
	  imass(3,3) = (2.0 *(eigenvalue[0][ind] - eigenvalue[3][ind] )) / Hartree /(k_max * k_max);
	  
	  
	  imass(2,1) = (eigenvalue[0][ind] - eigenvalue[4][ind] + (eigenvalue[1][ind] - eigenvalue[0][ind] )
	        + (eigenvalue[2][ind] - eigenvalue[0][ind]))  / Hartree /(k_max * k_max);



	  
	  imass(3,1) =(eigenvalue[0][ind] - eigenvalue[5][ind] + (eigenvalue[1][ind] - eigenvalue[0][ind] )
	         + (eigenvalue[3][ind] - eigenvalue[0][ind])) / Hartree /(k_max * k_max);
	 
	  imass(3,2) =(eigenvalue[0][ind] - eigenvalue[6][ind] + (eigenvalue[2][ind] - eigenvalue[0][ind] )
	         + (eigenvalue[3][ind] - eigenvalue[0][ind])) / Hartree /(k_max * k_max);
	 
	  double imass_DOS;
	  double temp1, temp2;
	  imass.invariants(&temp1, &temp2,&imass_DOS);

	 
	  
	  extremum.mass_DOS = std::pow(1.0/imass_DOS,1.0/3.0);
	  result.push_back(extremum);
	}
    }
  
 
  valence_band = result;

}
//----------------------------------------------------
void  WzDDsemiconductor::set_Ev(const double Ev)
{
  par.Ev = Ev;
}
//--------------------------------------------------
KPbulkHamiltonian::KPparams WzDDsemiconductor::calculate_6x6_kp_params (void )

{
 
  KPbulkHamiltonian::KPparams  result;



  //-----------------------------------------------
 // Valence band k.p parameters:

   
        

  result.L1 = 0.5 * (par.A5 + par.A4 + par.A2 - 1.0);
  result.L2 = 0.5 * (par.A1 - 1.0);
 
  result.M1 = 0.5 * (par.A4 + par.A2 - par.A5 - 1.0); 

  result.M2 = 0.5 * (par.A1 + par.A3 - 1.0);
        
  result.M3 = 0.5 * (par.A2 - 1.0);

  result.N1 = par.A5; 

  result.N2 = par.A6/sqrt(2.0);

  //------------------------------------------------------------------------------//
  // has to be checked!!!
  result.N1_yx = result.M1;  result.N2_yx = result.N1_yx;

  
  result.N1_xy = result.N1 - result.N1_yx; result.N2_xy = result.N1_xy;
  //-------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------//
 
  result.l1s =  par.D5/Hartree; result.l2s = par.D1/Hartree;

  result.m1s = (par.D4 + par.D2 - par.D5)/Hartree;  result.m2s = (par.D1 + par.D3)/Hartree;

  result.m3s = par.D2/Hartree;

  result.n1s = 2.0 * par.D5/Hartree; result.n2s = sqrt(2.0) * par.D6/Hartree;  
  //------------------------------------------------------------------------------//
  //------------------------------------------------------------------------------//
 
  

 

  //------------------------------------------------------------------------------//

  //spin-orbit energy
  result.d1 = par.delta_cr / Hartree ; //no crystal field splitting
  result.d2 = (par.delta_s/3) / Hartree;
  result.d3 = (par.delta_s/3) / Hartree;


  //valence band reference energy

  result.E_v = par.Ev/Hartree;
  
  
 
  


  return(result);
 
}
