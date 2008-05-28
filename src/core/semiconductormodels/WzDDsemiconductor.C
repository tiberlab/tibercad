#include "WzDDsemiconductor.h"
#include "WzSemiconductor.h"
#include "Constants.h"

using namespace std;
using namespace Constants;
 


 



//-------------------------------------------------------------/
void WzDDsemiconductor::do_calculate_conduction_band_extremum(void)
{
  vector<DDsemiconductor::band_extremum> result;
 
  double Ev_top;

  const WzSemiconductor::WzDDparameters& par = (dynamic_cast<WzSemiconductor*> (semiconductor))->get_parameters ();

  //---------------------------------
  //valence top reference energy
  //-------------------------------
  double d1 =  par.delta_cr;
  double d2 =  par.delta_s / 3.0;
  double d3 = d2;
  
   double E1 = d1 + d2;

   double E2 = (d1 - d2)/2.0 + sqrt( (d1-  d2)*( d1- d2) / 4.0 + 2.0 * d3 * d3 );

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
  band_ext.mass_DOS = pow(par.m_c_xx * par.m_c_xx * par.m_c_zz, 1.0/3.0) ;

  result.push_back(band_ext);

  conduction_band = result;
}
//----------------------------------------------------/

void  WzDDsemiconductor::do_calculate_valence_band_extremum(void)
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







