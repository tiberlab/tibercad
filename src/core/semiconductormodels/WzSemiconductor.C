// $Id$

#include "WzSemiconductor.h"
#include "Database.h"
#include "Material.h"

using namespace std;
 

const double WzSemiconductor::Hartree = 27.2113961;

 
//--------------------------------------------------//
WzSemiconductor::WzSemiconductor()
{


}
//--------------------------------------------------//
void WzSemiconductor::do_init()
{

  Semiconductor::do_init();
  
  ModelOptions & options = get_options ();
  {
    par.EgGamma      = get_parameter("Eg_G", par.EgGamma );
    par.Ev           = get_parameter("E_v", par.Ev );
    
    par.m_c_zz       = options.get_option("m_c_zz", par.m_c_zz);
    par.m_c_xx       = options.get_option("m_c_xx", par.m_c_xx);
  
    par.A1           = options.get_option("A1", par.A1);
    par.A2           = options.get_option("A2", par.A2);
    par.A3           = options.get_option("A3", par.A3);
    par.A4           = options.get_option("A4", par.A4);
    par.A5           = options.get_option("A5", par.A5);
    par.A6           = options.get_option("A6", par.A6); 
  
    par.a_x          = options.get_option("a_x", par.a_x);
    par.a_z          = options.get_option("a_z", par.a_z);
    
    par.D1           = options.get_option("D1", par.D1 );
    par.D2           = options.get_option("D2",par.D2 );
    par.D3           = options.get_option("D3", par.D3 );
    par.D4           = options.get_option("D4",  par.D4);
    par.D5           = options.get_option("D5", par.D5);
    par.D6           = options.get_option("D6", par.D6 );
    par.delta_s      = options.get_option("delta_s",  par.delta_s);
    par.delta_cr     = options.get_option("delta_cr",  par.delta_cr);

  
    par.Ep_1 = options.get_option("Ep_1", par.Ep_1);
    par.Ep_2 = options.get_option("Ep_2", par.Ep_2);

    par.varshni_alpha_G = options.get_option("varshni_alpha_G", par.varshni_alpha_G );
    par.varshni_beta_G = options.get_option("varshni_beta_G",  par.varshni_beta_G);
   


  }

  {
    bow.EgGamma      = get_parameter("bow_Eg_G",bow.EgGamma );
    bow.Ev           = options.get_option("bow_E_v", bow.Ev );
    
    bow.m_c_zz       = options.get_option("bow_m_c_zz", bow.m_c_zz);
    bow.m_c_xx       = options.get_option("bow_m_c_xx", bow.m_c_xx);
  
    bow.A1           = options.get_option("bow_A1", bow.A1);
    bow.A2           = options.get_option("bow_A2", bow.A2);
    bow.A3           = options.get_option("bow_A3", bow.A3);
    bow.A4           = options.get_option("bow_A4", bow.A4);
    bow.A5           = options.get_option("bow_A5", bow.A5);
    bow.A6           = options.get_option("bow_A6", bow.A6); 
  
    bow.a_x          = options.get_option("bow_a_x", bow.a_x);
    bow.a_z          = options.get_option("bow_a_z", bow.a_z);
    
    bow.D1           = options.get_option("bow_D1", bow.D1 );
    bow.D2           = options.get_option("bow_D2",bow.D2 );
    bow.D3           = options.get_option("bow_D3", bow.D3 );
    bow.D4           = options.get_option("bow_D4",  bow.D4);
    bow.D5           = options.get_option("bow_D5", bow.D5);
    bow.D6           = options.get_option("bow_D6", bow.D6 );
    bow.delta_s      = options.get_option("bow_delta_s",  bow.delta_s);
    bow.delta_cr     = options.get_option("bow_delta_cr",  bow.delta_cr);

  
    bow.Ep_1 = options.get_option("bow_Ep_1", bow.Ep_1);
    bow.Ep_2 = options.get_option("bow_Ep_2", bow.Ep_2);
  }
 


   {
    //here zero temperature and work parameters coinside
    par_initial = par;
   }


}

//--------------------------------------------------//
void  WzSemiconductor::read_database_alloy(void)
{

  Database& db = get_database();

  db.set_section("bandgap");
  bow.EgGamma = db.get("bow_Eg_G", 0.0);


  db.set_section("kdotp");
  bow.Ev = db.get("bow_E_v", 0.0);

  bow.m_c_zz = db.get("bow_m_c_zz", 0.0);
  bow.m_c_xx = db.get("bow_m_c_xx", 0.0);
  
  bow.A1 = db.get("bow_A1", 0.0);
  bow.A2 = db.get("bow_A2", 0.0);
  bow.A3 = db.get("bow_A3", 0.0);
  bow.A4 = db.get("bow_A4", 0.0);
  bow.A5 = db.get("bow_A5", 0.0);
  bow.A6 = db.get("bow_A6", 0.0); 

  bow.delta_s = db.get("bow_delta_s", 0.0);
  bow.delta_cr = db.get("bow_delta_cr", 0.0);
  
  bow.Ep_1 = db.get("bow_Ep_1", 0.0);
  bow.Ep_2 = db.get("bow_Ep_2", 0.0);
  
  db.set_section("deformation_potentials");
  bow.a_x = db.get("bow_a_x", 0.0);
  bow.a_z = db.get("bow_a_z", 0.0);
  
  bow.D1 = db.get("bow_D1", 0.0);
  bow.D2 = db.get("bow_D2", 0.0);
  bow.D3 = db.get("bow_D3", 0.0);
  bow.D4 = db.get("bow_D4", 0.0);
  bow.D5 = db.get("bow_D5", 0.0);
  bow.D6 = db.get("bow_D6", 0.0);

}


//--------------------------------------------------

void WzSemiconductor::read_database( )
{
  Database& db = get_database();

  db.set_section("bandgap");
  par.EgGamma = db.get("Eg_G", 3.51);

  par.varshni_alpha_G = db.get("varshni_alpha_G", 0.0);
  par.varshni_beta_G  = db.get("varshni_beta_G", 0.0);


  db.set_section("kdotp");
  par.Ev = db.get("E_v", -0.726);

  par.m_c_zz = db.get("m_c_zz", 0.20);
  par.m_c_xx = db.get("m_c_xx", 0.20);
  
  par.A1 = db.get("A1", -7.21);
  par.A2 = db.get("A2", -0.44);
  par.A3 = db.get("A3", 6.68);
  par.A4 = db.get("A4", -3.46);
  par.A5 = db.get("A5", -3.40);
  par.A6 = db.get("A6", -4.90); 

  par.delta_s = db.get("delta_s", 0.017);
  par.delta_cr = db.get("delta_cr", 0.010);
  
  par.Ep_1 = db.get("Ep_1", 14.0);
  par.Ep_2 = db.get("Ep_2", 14.0);
  

  db.set_section("deformation_potentials");
  par.a_x = db.get("a_x", -4.9);
  par.a_z = db.get("a_z", -11.3);
  
  par.D1 = db.get("D1", -3.7);
  par.D2 = db.get("D2", 4.5);
  par.D3 = db.get("D3", 8.2);
  par.D4 = db.get("D4", -4.1);
  par.D5 = db.get("D5", -4.0);
  par.D6 = db.get("D6", -5.5);



  {
    //here zero temperature and work parameters coinside
    par_initial = par;
  }

 

}


//--------------------------------------------------------------------//
void WzSemiconductor::do_do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
 

  const WzSemiconductor* modA = dynamic_cast<const WzSemiconductor*> (comp_A);
  const WzSemiconductor* modB = dynamic_cast<const WzSemiconductor*> (comp_B);
  
  par.EgGamma      = alloy(modA->par.EgGamma, modB->par.EgGamma, xa, bow.EgGamma);
  par.Ev           = alloy(modA->par.Ev, modB->par.Ev, xa, bow.Ev);
    
  par.m_c_zz       = alloy(modA->par.m_c_zz, modB->par.m_c_zz, xa, bow.m_c_zz);
  par.m_c_xx       = alloy(modA->par.m_c_xx, modB->par.m_c_xx, xa, bow.m_c_xx);
  
  par.A1           = alloy( modA->par.A1, modB->par.A1, xa, bow.A1);
  par.A2           = alloy( modA->par.A2, modB->par.A2, xa, bow.A2);
  par.A3           = alloy( modA->par.A3, modB->par.A3, xa, bow.A3);
  par.A4           = alloy( modA->par.A4, modB->par.A4, xa, bow.A4);
  par.A5           = alloy( modA->par.A5, modB->par.A5, xa, bow.A5);
  par.A6           = alloy( modA->par.A6, modB->par.A6, xa, bow.A6);
  
  par.a_x          = alloy(modA->par.a_x, modB->par.a_x, xa, bow.a_x);
  par.a_z          = alloy(modA->par.a_z, modB->par.a_z, xa, bow.a_z);
    
  par.D1           = alloy( modA->par.D1, modB->par.D1, xa, bow.D1);
  par.D2           = alloy( modA->par.D2, modB->par.D2, xa, bow.D2);
  par.D3           = alloy( modA->par.D3, modB->par.D3, xa, bow.D3);
  par.D4           = alloy( modA->par.D4, modB->par.D4, xa, bow.D4);
  par.D5           = alloy( modA->par.D5, modB->par.D5, xa, bow.D5);
  par.D6           = alloy( modA->par.D6, modB->par.D6, xa, bow.D6);

  par.delta_s      = alloy(modA->par.delta_s, modB->par.delta_s, xa, bow.delta_s);
  par.delta_cr     = alloy(modA->par.delta_cr, modB->par.delta_cr, xa, bow.delta_cr);

  
  par.Ep_1 = alloy(modA->par.Ep_1, modB->par.Ep_1, xa, bow.Ep_1);
  par.Ep_2 = alloy(modA->par.Ep_2, modB->par.Ep_2, xa, bow.Ep_2);


  {
    //here zero temperature and work parameters coinside
    par_initial = par;
  }



}
//--------------------------------------------------
KPparams WzSemiconductor::do_calculate_8x8_kp_params (void )
{
  //we start from 6x6 kp parameters  
  KPparams  result = calculate_6x6_kp_params();

  
  

  //------------------------------------------------------------------
  //CONDUCTION BAND
  //we renormalize conduction band quadratic part to free electron mass
  result.s1 = 1.0;
  result.s2 = 1.0;


  //---------------------------------
  //valence top reference energy
  //-------------------------------

  double Ev_top;


  {
    double d1 =  par.delta_cr;
    double d2 =  par.delta_s / 3.0;
    double d3 = d2;
  
    double E1 = d1 + d2;

    double E2 = (d1 - d2)/2.0 + sqrt( (d1 -  d2)*( d1 - d2) / 4.0 + 2.0 * d3 * d3 );

    if (E1 > E2)
      Ev_top = par.Ev + E1;
    else
      Ev_top = par.Ev + E2;

 


   
  }

  result.E_c =  (Ev_top + par.EgGamma)/Hartree;



  //-----------------------------------------------------------------
  //to check!
  double Ep1,Ep2; //Ep = P^2 * 2.0 (atomic units);

  double s = 1.0;
  
  Ep1 = s*(1.0/par.m_c_zz - 1.0)*
    par.EgGamma * ( (par.EgGamma +  par.delta_s )/(par.EgGamma + 2.0/3.0 * par.delta_s ) ) / Hartree;

  Ep2 = s*(1.0/par.m_c_xx - 1.0)*
    par.EgGamma * ( (par.EgGamma +  par.delta_s )/(par.EgGamma + 2.0/3.0 * par.delta_s ) ) / Hartree;
  

  result.P1 = std::sqrt(0.5 * Ep1);
  result.P2 = std::sqrt(0.5 * Ep2);

 


  //-----------------------------------------------------------------
  //to check !
  //rescale L and N
  double t1 =   0.5*Ep1/( (par.EgGamma +  par.delta_s/3.0)/Hartree );
  double t2 =   0.5*Ep2/( (par.EgGamma +  par.delta_s/3.0)/Hartree );

  result.L1 += t1;
  result.L2 += t2;

  result.N1 += t1;
  result.N2 += t2;
  

  result.N1_yx = result.M1;  result.N2_yx = result.N1_yx;
  result.N1_xy = result.N1 - result.N1_yx; result.N2_xy = result.N1_xy;

  return(result);



}

//--------------------------------------------------
KPparams WzSemiconductor::do_calculate_6x6_kp_params (void )

{
 
  KPparams  result;



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
  result.N1_yx = result.M1;  result.N2_yx = result.M2;

  
  result.N1_xy = result.N1 - result.N1_yx; result.N2_xy = result.N2 - result.N2_yx;
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
  result.P1 = std::sqrt(par.Ep_1 * 0.5 / Hartree);
  result.P2 = std::sqrt(par.Ep_2 * 0.5 / Hartree);


  return(result);
 
}


//---------------------------
void WzSemiconductor::apply_temperature(void)
{

  par = par_initial;


  
  if (get_material()->is_alloy())
  {
    const WzSemiconductor::WzDDparameters& parA_zero =  (dynamic_cast<const WzSemiconductor*> (modelA))->get_initial_parameters();

    const WzSemiconductor::WzDDparameters& parB_zero =  (dynamic_cast<const WzSemiconductor*> (modelB))->get_initial_parameters();

     double EgGammaA = parA_zero.EgGamma 
      - parA_zero.varshni_alpha_G * _temperature*_temperature/(_temperature + parA_zero.varshni_beta_G);
    double EgGammaB = parB_zero.EgGamma 
      - parB_zero.varshni_alpha_G * _temperature*_temperature/(_temperature + parB_zero.varshni_beta_G);

    par.EgGamma = alloy( EgGammaA, EgGammaB, _xa, bow.EgGamma);

  }
  else
  {
    par.EgGamma = par_initial.EgGamma - par.varshni_alpha_G * _temperature*_temperature/(_temperature + par.varshni_beta_G) ;
  }
}
