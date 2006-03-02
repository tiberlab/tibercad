using namespace std;
#include "ZbDDsemiconductor.h"
//---------------------------------------------//

ZbDDsemiconductor::ZbDDsemiconductor(void): DDsemiconductor()
{
  EgGamma = 0;
  EgL = 0;
  EgX = 0;
  gamma1 = 0;
  gamma2 = 0;
  gamma3 = 0; 
  m_t_L = 0;
  m_l_L = 0;    
  m_t_X = 0;
  m_l_X = 0;  
  m_G = 0;             
  a_c = 0;
  a_v = 0;
  b = 0;
  d = 0; 
  def_vol_X = 0; 
  def_uniax_X = 0; //volume and uniax deformation potential for X point
  def_vol_L = 0;
  def_uniax_L = 0; //volume and uniax deformation potential for L point
  delta = 0;
} 

//--------------------------------------------//
ZbDDsemiconductor::ZbDDsemiconductor(const ZbDDparameters& params): DDsemiconductor()
{
  EgGamma = params.EgGamma;
  EgL = params.EgL;
  EgX = params.EgX ;
  gamma1 = params.gamma1;
  gamma2 = params.gamma2;
  gamma3 = params.gamma3; 
  m_t_L = params.m_t_L;
  m_l_L = params.m_l_L;    
  m_t_X = params.m_t_X;
  m_l_X = params.m_l_X;  
  m_G = params.m_G;             
  a_c = params.a_c;
  a_v = params.a_v;
  b = params.b;
  d = params.d; 
  def_vol_X =  params.def_vol_X ; 
  def_uniax_X =  params.def_uniax_X; 
  def_vol_L =  params.def_vol_L;
  def_uniax_L =  params.def_uniax_L;
}

//---------------------------------------------//

void ZbDDsemiconductor::set_Eg(double EgGamma_1, double EgL_1, double EgX_1)
{
  EgGamma = EgGamma_1;
  EgL = EgL_1;
  EgX = EgX_1;

}

//--------------------------------------------//

void ZbDDsemiconductor::set_mass_Gamma(double m)
{
  m_G = m;
}

//--------------------------------------------//

void ZbDDsemiconductor::set_masses_L(double m_t, double m_l)
{
  m_t_L = m_t;
  m_l_L = m_l;
}

//-------------------------------------------//
void ZbDDsemiconductor::set_masses_X(double m_t, double m_l)
{
  m_t_X = m_t;
  m_l_X = m_l;
}

//-------------------------------------------//

void ZbDDsemiconductor::set_6x6kp_params(double gamma1_, double gamma2_, double gamma3_, double delta_)
{
  gamma1 = gamma1_;
  gamma2 = gamma2_;
  gamma3 = gamma3_;
  delta  = delta_;
}

//--------------------------------------------//
void ZbDDsemiconductor::set_deformation_parameters(double a_c_, double a_v_, double b_, double d_)
{
  a_c = a_c_;
  a_v = a_v_;
  b = b_;
  d = d_;
}
//--------------------------------------------------------//
void ZbDDsemiconductor::calculate_conduction_band_extremum(void)
{
  vector<DDsemiconductor::band_extremum>   result;
  band_extremum  band_ex;
  //there are 3 types of extrema: Gamma, X and L
  //---------------------------------------------------------
  //Gamma minima
  //---------------------------------------------------------

  double Ec_G = Ev + EgGamma;
 

  if (strained)  Ec_G  += a_c * trace(strain);
    
 
  band_ex.degeneracy = 2   ; //spin
  band_ex.energy     = Ec_G;
  band_ex.mass_DOS   = m_G ; //mass is isotropic

  result.push_back(band_ex);
  //----------------------------------------------------------------
  //L minima
  //----------------------------------------------------------------
  double Ec_L = Ev + EgL;
  if (strained)
    {
      
      //Hydrostatic strain part--------------------------------
      Ec_L  += def_vol_L * trace(strain);
      //Uniaxial strain part----------------------------------
      Tensor1 k;
      vector<double> dE_uniax(4);
      k(1) = 1 ; k(2) =  1; k(3) =  1;
      dE_uniax[0] = def_uniax_L * (k * (strain * k)) * (1.0/3.0);

      k(1) = -1; k(2) =  1; k(3) =  1;
      dE_uniax[1] = def_uniax_L * (k * (strain * k)) * (1.0/3.0);

     

      k(1) =  1; k(2) = -1; k(3) =  1;
      dE_uniax[2] = def_uniax_L * (k * (strain * k)) * (1.0/3.0);

      

      k(1) =  1; k(2) =  1; k(3) = -1;
      dE_uniax[3] = def_uniax_L * (k * (strain * k)) * (1.0/3.0);

     
      
      double mass_DOS = pow(m_t_L * m_t_L *  m_l_L, 1.0/3.0 );
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
      band_ex.degeneracy = 16   ; /* spin degeneracy and  8 equivalent minima (
                                  [1,1, 1],[-1,-1,-1], 
                                  [-1,1,1],[1,-1,-1] 
				  [1,-1,1],[-1,1,-1]
                                  [1,1,-1],[-1,-1,1]
				  */

      band_ex.energy     = Ec_L;
      band_ex.mass_DOS = pow(m_t_L * m_t_L *  m_l_L, 1.0/3.0 );

      result.push_back(band_ex);
    }

  //------------------------------------------------------------------
  //X-minima
  //------------------------------------------------------------------
  double Ec_X = Ev + EgX;
  if (strained)
    {
      //Hydrostatic strain part--------------------------------
      Ec_X  += def_vol_X * trace(strain);
      //-------------------------------------------------------
      //Uniaxial strain part----------------------------------
      Tensor1 k;
      vector<double> dE_uniax(3);
      k(1) = 1 ; k(2) =  0; k(3) =  0;
      dE_uniax[0] = def_uniax_X * (k * (strain * k));

      k(1) = 0 ; k(2) =  1; k(3) =  0;
      dE_uniax[1] = def_uniax_X * (k * (strain * k));

      k(1) = 0 ; k(2) =  0; k(3) =  1;
      dE_uniax[2] = def_uniax_X * (k * (strain * k));
      
      double mass_DOS = pow(m_t_X * m_t_X *  m_l_X, 1.0/3.0 );
      
      for (short i = 0; i <3; i++)
	{  
	  band_ex.degeneracy = 4; //spin and k-> -k
	  band_ex.mass_DOS = mass_DOS;
	  band_ex.energy = Ec_X + dE_uniax[i];
	  result.push_back(band_ex);
	}

      
    }
  else
    {
      band_ex.energy     = Ec_X;
      band_ex.degeneracy = 12   ; // spin degeneracy and 6 equivalent minima
      band_ex.mass_DOS = pow(m_t_X * m_t_X *  m_l_X, 1.0/3.0 );
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

//--------------------------------------------------------//

void ZbDDsemiconductor::calculate_valence_band_energy_extremum(void)
{

}
//---------------------------------------------------------//
