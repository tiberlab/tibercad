#include "ZbDDsemiconductor.h"
#include "ZbSemiconductor.h"
#include "getpot.h"
#include "Alloy.h"
#include "Database.h"

using namespace std;


//---------------------------------------------//


void ZbDDsemiconductor::calculate_conduction_band_extremum(void)

{
   vector<DDsemiconductor::band_extremum>   result;

   const ZbSemiconductor::ZbDDparameters& par = (dynamic_cast<ZbSemiconductor*> (semiconductor))->get_parameters ();

   band_extremum  band_ex;
   //there are 3 types of extrema: Gamma, X and L
   //---------------------------------------------------------
   //Gamma minima
   //---------------------------------------------------------

  

   double Ev_top = par.Ev + ((1.0/3.0) * par.delta );

   double Ec_G = Ev_top + par.EgGamma;
 

   if (strained)  Ec_G  += par.a_c * trace(strain);
    
 
   band_ex.degeneracy = 2   ; //spin
   band_ex.energy     = Ec_G;
   band_ex.mass_DOS   = par.m_G ; //mass is isotropic

   result.push_back(band_ex);
   //----------------------------------------------------------------
   //L minima
   //----------------------------------------------------------------
   double Ec_L = Ev_top + par.EgL;
   if (strained)
     {
      
       //Hydrostatic strain part--------------------------------
       Ec_L  += (par.def_vol_L - (1.0/3.0)*par.def_uniax_L ) * trace(strain);
       //Uniaxial strain part----------------------------------
       Tensor1 k;
       vector<double> dE_uniax(4);
       k(1) = 1 ; k(2) =  1; k(3) =  1;
       dE_uniax[0] = par.def_uniax_L * (k * (strain * k)) * (1.0/3.0);

       k(1) = -1; k(2) =  1; k(3) =  1;
       dE_uniax[1] = par.def_uniax_L * (k * (strain * k)) * (1.0/3.0);

     

       k(1) =  1; k(2) = -1; k(3) =  1;
       dE_uniax[2] = par.def_uniax_L * (k * (strain * k)) * (1.0/3.0);

      

       k(1) =  1; k(2) =  1; k(3) = -1;
       dE_uniax[3] = par.def_uniax_L * (k * (strain * k)) * (1.0/3.0);

     
      
       double mass_DOS = pow(par.m_t_L * par.m_t_L *  par.m_l_L, 1.0/3.0 );
       for (short i = 0; i <4; i++)
	 {
	   band_ex.degeneracy = 4; //spin and k-> -k
	   band_ex.mass_DOS = mass_DOS;
	  
	   band_ex.energy = Ec_L + dE_uniax[i];
	   result.push_back(band_ex);
	 }

     }
   else
     {
       band_ex.degeneracy = 16   ; 
       /* spin degeneracy and  8 equivalent minima (
	  [1,1, 1],[-1,-1,-1], 
	  [-1,1,1],[1,-1,-1] 
	  [1,-1,1],[-1,1,-1]
	  [1,1,-1],[-1,-1,1]
       */

       band_ex.energy     = Ec_L;
       band_ex.mass_DOS = pow(par.m_t_L * par.m_t_L *  par.m_l_L, 1.0/3.0 );

       result.push_back(band_ex);
     }
 
   //------------------------------------------------------------------
   //X-minima
   //------------------------------------------------------------------
   double Ec_X = Ev_top + par.EgX;
   if (strained)
   {
     //Hydrostatic strain part--------------------------------
     Ec_X  += (par.def_vol_X - (1.0/3.0)*par.def_uniax_X) * trace(strain);
     //-------------------------------------------------------
     //Uniaxial strain part----------------------------------
     Tensor1 k;
     vector<double> dE_uniax(3);
     k(1) = 1 ; k(2) =  0; k(3) =  0;
     dE_uniax[0] = par.def_uniax_X * (k * (strain * k));
     
     k(1) = 0 ; k(2) =  1; k(3) =  0;
     dE_uniax[1] = par.def_uniax_X * (k * (strain * k));
     
     k(1) = 0 ; k(2) =  0; k(3) =  1;
     dE_uniax[2] = par.def_uniax_X * (k * (strain * k));
     
     double mass_DOS = pow(par.m_t_X * par.m_t_X *  par.m_l_X, 1.0/3.0 );
     
     for (short i = 0; i <3; i++)
     {  
       band_ex.degeneracy = 4; //spin and k->-k
       band_ex.mass_DOS = mass_DOS;
       band_ex.energy = Ec_X + dE_uniax[i];
       result.push_back(band_ex);
     }

      
   }
   else
   {
     band_ex.energy     = Ec_X;
     band_ex.degeneracy = 12   ; // spin degeneracy and 6 equivalent minima
     band_ex.mass_DOS = pow(par.m_t_X * par.m_t_X *  par.m_l_X, 1.0/3.0 );
     result.push_back(band_ex);
   }
   
   //------------------------------------------
   //let us find the minimal energy
   short N = result.size();
   double Emin = result[0].energy;
   for (short i = 1; i < N; i++)
   {
     if (Emin > result[i].energy) Emin = result[i].energy;
   }
   //-----------------------------------------
   
   vector<DDsemiconductor::band_extremum> result_filtered;
   for (short i = 0; i < N; i++)
   {
     if ( result[i].energy - Emin <= energy_cutoff)
     {
       result_filtered.push_back(result[i]);
     }
   }

   //-----------------------------------------
  

   conduction_band = result_filtered;

}





//=================================================================================//


void ZbDDsemiconductor::calculate_valence_band_extremum(void)
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

  for (short ind = 0; ind < 3; ind++)
  {
    if (eigenvalue[0][ind*2] + energy_cutoff > eigenvalue[0][5]) 
    {
      extremum.degeneracy = 2;
      extremum.energy  =eigenvalue[0][ind*2] ;
	  
      imass(1,1) = (2.0 *(eigenvalue[0][ind*2] - eigenvalue[1][ind*2] )) / Hartree /(k_max * k_max);

	 
      imass(2,2) = (2.0 *(eigenvalue[0][ind*2] - eigenvalue[2][ind*2] )) / Hartree /(k_max * k_max);
	 
      imass(3,3) = (2.0 *(eigenvalue[0][ind*2] - eigenvalue[3][ind*2] )) / Hartree /(k_max * k_max);
	  
	  
      imass(2,1) = (eigenvalue[0][ind*2] - eigenvalue[4][ind*2] + (eigenvalue[1][ind*2] - eigenvalue[0][ind*2] )
		    + (eigenvalue[2][ind*2] - eigenvalue[0][ind*2]))  / Hartree /(k_max * k_max);



	  
      imass(3,1) =(eigenvalue[0][ind*2] - eigenvalue[5][ind*2] + (eigenvalue[1][ind*2] - eigenvalue[0][ind*2] )
		   + (eigenvalue[3][ind*2] - eigenvalue[0][ind*2])) / Hartree /(k_max * k_max);
	 
      imass(3,2) =(eigenvalue[0][ind*2] - eigenvalue[6][ind*2] + (eigenvalue[2][ind*2] - eigenvalue[0][ind*2] )
		   + (eigenvalue[3][ind*2] - eigenvalue[0][ind*2])) / Hartree /(k_max * k_max);
	 
      double imass_DOS;
      double temp1, temp2;
      

      imass.invariants(&temp1, &temp2,&imass_DOS);

	  
      extremum.mass_DOS = std::pow(1.0/imass_DOS,1.0/3.0);
      result.push_back(extremum);
    }
  }
  
 
  valence_band = result;
}

//---------------------------------------------------------//


//============================================================================================//
