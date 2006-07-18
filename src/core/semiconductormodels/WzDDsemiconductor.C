using namespace std;
#include "WzDDsemiconductor.h"

#include "Alloy.h"
#include "getpot.h"

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
  par.Ep_1 = 0;
  par.Ep_2 = 0;
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

void WzDDsemiconductor::read_database(const Dummy&)
{

  GetPot data(_filename);

  const std::string structure = data("structure", "wz");

  assert(structure == "wz");


  par.EgGamma = data("Eg_G", 3.51);
  par.Ev = data("E_v", -0.726);

  par.m_c_zz = data("m_c_zz", 0.20);
  par.m_c_xx = data("m_c_xx", 0.20);
  
  par.A1 = data("A1", -7.21);
  par.A2 = data("A2", -0.44);
  par.A3 = data("A3", 6.68);
  par.A4 = data("A4", -3.46);
  par.A5 = data("A5", -3.40);
  par.A6 = data("A6", -4.90); 
  
  par.a_x = data("a_x", -4.9);
  par.a_z = data("a_z", -11.3);
  
  par.D1 = data("D1", -3.7);
  par.D2 = data("D2", 4.5);
  par.D3 = data("D3", 8.2);
  par.D4 = data("D4", -4.1);
  par.D5 = data("D5", -4.0);
  par.D6 = data("D6", -5.5);
  par.delta_s = data("delta_s", 0.017);
  par.delta_cr = data("delta_cr", 0.010);

  
  par.Ep_1 = data("Ep_1", 14.0);
  par.Ep_2 = data("Ep_2", 14.0);


}


//--------------------------------------------------
void WzDDsemiconductor::build_alloy(const std::string& component2,
		 const std::string& bowing_params, double content)
{

  GetPot data(component2);
  GetPot bowing(bowing_params);

  const std::string structure = data("structure", "wz");

  assert(structure == "wz");

  double (*alloy)(double, double, double, double) =
    Alloy::calculate_VCA_parameter;
 

  par.EgGamma = alloy(data("Eg_G", 3.51), par.EgGamma, content,
			 bowing("Eg_G", 0.0));
  par.Ev = alloy(data("E_v", -0.726), par.Ev, content,
		    bowing("E_v", 0.0));

  par.m_c_zz = alloy(data("m_c_zz", 0.20), par.m_c_zz, content,
			bowing("m_c_zz", 0.0));
  par.m_c_xx = alloy(data("m_c_xx", 0.20), par.m_c_xx, content,
			bowing("m_c_xx", 0.0));
  
  par.A1 = alloy(data("A1", -7.21), par.A1, content,
		    bowing("A1", 0.0));
  par.A2 = alloy(data("A2", -0.44), par.A2, content,
		    bowing("A2", 0.0));
  par.A3 = alloy(data("A3", 6.68), par.A3, content,
		    bowing("A3", 0.0));
  par.A4 = alloy(data("A4", -3.46), par.A4,content, 
		    bowing("A4", 0.0)); 
  par.A5 = alloy(data("A5", -3.40), par.A5, content,
		    bowing("A5", 0.0));
  par.A6 = alloy(data("A6", -4.90), par.A6,content, 
		    bowing("A6", 0.0));
  
  par.a_x = alloy(data("a_x", -4.9), par.a_x, content,
		     bowing("a_x", 0.0));
  par.a_z = alloy(data("a_z", -11.3), par.a_z,content, 
		     bowing("a_z", 0.0));
  
  par.D1 = alloy(data("D1", -3.7), par.D1, content,
		    bowing("D1", 0.0));
  par.D2 = alloy(data("D2", 4.5), par.D2,content, 
		    bowing("D2", 0.0));
  par.D3 = alloy(data("D3", 8.2), par.D3, content,
		    bowing("D3", 0.0));
  par.D4 = alloy(data("D4", -4.1), par.D4, content,
		    bowing("D4", 0.0));
  par.D5 = alloy(data("D5", -4.0), par.D5, content,
		    bowing("D5", 0.0));
  par.D6 = alloy(data("D6", -5.5), par.D6, content,
		    bowing("D6", 0.0));
  par.delta_s = alloy(data("delta_s", 0.017), par.delta_s,
			 content, bowing("delta_s", 0.0));
  par.delta_cr = alloy(data("delta_cr", 0.010), par.delta_cr,
			  content, bowing("delta_cr", 0.0));

  
  par.Ep_1 = alloy(data("Ep_1", 14.0), par.delta_cr,
			  content, bowing("Ep_1", 0.0));

  par.Ep_2 = alloy(data("Ep_2", 14.0), par.delta_cr,
			  content, bowing("Ep_2", 0.0)); 
  
}
//--------------------------------------------------
KPbulkHamiltonian::KPparams WzDDsemiconductor::calculate_8x8_kp_params (void )
{


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
 
  result.l1s =  (par.D5 +  par.D4 + par.D2)/Hartree; result.l2s = par.D1/Hartree;

  result.m1s = (par.D4 + par.D2 - par.D5)/Hartree;  result.m2s = (par.D1 + par.D3)/Hartree;

  result.m3s = par.D2/Hartree;

  result.n1s = 2.0 * par.D5/Hartree; result.n2s = sqrt(2.0) * par.D6/Hartree;  
  //------------------------------------------------------------------------------//
  //------------------------------------------------------------------------------//
  // conduction band strain
  result.axs = par.a_x / Hartree;
  result.azs = par.a_z / Hartree;
  
 

  //------------------------------------------------------------------------------//

  //spin-orbit energy
  result.d1 = par.delta_cr / Hartree ; //no crystal field splitting
  result.d2 = (par.delta_s/3) / Hartree;
  result.d3 = (par.delta_s/3) / Hartree;


  //valence band reference energy

  result.E_v = par.Ev/Hartree;
  

 
  
  
  //conduction-valence band coupling. may be needed only for optics
  result.P1 = std::sqrt(par.Ep_1 * 2.0 / Hartree);
  result.P2 = std::sqrt(par.Ep_2 * 2.0 / Hartree);


  return(result);
 
}
