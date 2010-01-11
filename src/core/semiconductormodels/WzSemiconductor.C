// $Id$

#include "WzSemiconductor.h"
#include "Database.h"
#include "Material.h"
#include "Constants.h"

using namespace std;



//--------------------------------------------------//
WzSemiconductor::WzSemiconductor(const ModelOptions& options)
 : Semiconductor(options)
{


}
//--------------------------------------------------//
void WzSemiconductor::do_init()
{

  Semiconductor::do_init();

  ModelOptions & options = get_options ();
  {
    get_parameter("Eg_G", par.EgGamma );
    get_parameter("E_v", par.Ev );

    par.m_c_zz       = get_option("m_c_zz", par.m_c_zz);
    par.m_c_xx       = get_option("m_c_xx", par.m_c_xx);

    par.A1           = get_option("A1", par.A1);
    par.A2           = get_option("A2", par.A2);
    par.A3           = get_option("A3", par.A3);
    par.A4           = get_option("A4", par.A4);
    par.A5           = get_option("A5", par.A5);
    par.A6           = get_option("A6", par.A6);

    par.a_x          = get_option("a_x", par.a_x);
    par.a_z          = get_option("a_z", par.a_z);

    par.D1           = get_option("D1", par.D1 );
    par.D2           = get_option("D2",par.D2 );
    par.D3           = get_option("D3", par.D3 );
    par.D4           = get_option("D4",  par.D4);
    par.D5           = get_option("D5", par.D5);
    par.D6           = get_option("D6", par.D6 );
    par.delta_s      = get_option("delta_s",  par.delta_s);
    par.delta_cr     = get_option("delta_cr",  par.delta_cr);


    par.Ep_1 = get_option("Ep_1", par.Ep_1);
    par.Ep_2 = get_option("Ep_2", par.Ep_2);

    par.varshni_alpha_G = get_option("varshni_alpha_G", par.varshni_alpha_G );
    par.varshni_beta_G = get_option("varshni_beta_G",  par.varshni_beta_G);



  }
}


//--------------------------------------------------//
void  WzSemiconductor::read_database_alloy(void)
{

  Database& db = get_database();

  db.set_section("bandgap");
  bow.EgGamma = db.get("bow_Eg_G", 0.0);

}


//--------------------------------------------------

void WzSemiconductor::read_database( )
{
  Database& db = get_database();

  db.set_section("bandgap");
  par.EgGamma = db.get("Eg_G", 3.51);

  par.varshni_alpha_G = db.get("varshni_alpha_G", 0.0);
  par.varshni_beta_G  = db.get("varshni_beta_G", 0.0);


  db.set_section("valenceband");
  par.Ev = db.get("E_v", -0.726);

  db.set_section("conductionband");
  par.m_c_zz = db.get("m_c_zz", 0.20);
  par.m_c_xx = db.get("m_c_xx", 0.20);

  db.set_section("kdotp");
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

  result.E_c =  (par.Ev + par.EgGamma)/Constants::Hartree;



  //-----------------------------------------------------------------
  //to check!
  double Ep1,Ep2; //Ep = P^2 * 2.0 (atomic units);

  double s = 1.0;

  Ep1 = s*(1.0/par.m_c_zz - 1.0)*
    par.EgGamma * ( (par.EgGamma +  par.delta_s )/(par.EgGamma + 2.0/3.0 * par.delta_s ) ) / Constants::Hartree;

  Ep2 = s*(1.0/par.m_c_xx - 1.0)*
    par.EgGamma * ( (par.EgGamma +  par.delta_s )/(par.EgGamma + 2.0/3.0 * par.delta_s ) ) / Constants::Hartree;


  result.P1 = std::sqrt(0.5 * Ep1);
  result.P2 = std::sqrt(0.5 * Ep2);




  //-----------------------------------------------------------------
  //to check !
  //rescale L and N
  double t1 =   0.5*Ep1/( (par.EgGamma +  par.delta_s/3.0)/Constants::Hartree );
  double t2 =   0.5*Ep2/( (par.EgGamma +  par.delta_s/3.0)/Constants::Hartree );

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

  result.l1s =  (par.D5 +  par.D4 + par.D2)/Constants::Hartree; result.l2s = par.D1/Constants::Hartree;

  result.m1s = (par.D4 + par.D2 - par.D5)/Constants::Hartree;  result.m2s = (par.D1 + par.D3)/Constants::Hartree;

  result.m3s = par.D2/Constants::Hartree;

  result.n1s = 2.0 * par.D5/Constants::Hartree; result.n2s = sqrt(2.0) * par.D6/Constants::Hartree;


  //------------------------------------------------------------------------------//
  // conduction band strain
  result.axs = par.a_x / Constants::Hartree;
  result.azs = par.a_z / Constants::Hartree;



  //------------------------------------------------------------------------------//
  //spin-orbit energy
  result.d1 = par.delta_cr / Constants::Hartree ;
  result.d2 = (par.delta_s/3) / Constants::Hartree;
  result.d3 = (par.delta_s/3) / Constants::Hartree;


  //------------------------------------------------------------------------------//
  //valence band reference energy
  {
    double d1 = result.d1;
    double d2 = result.d2;
    double d3 = result.d3;

    double E1 = d1 + d2;
    double E2 = (d1 - d2)/2.0
      + sqrt( (d1 -  d2)*( d1 - d2) / 4.0 + 2.0 * d3 * d3 );

    if (E1 > E2)
      result.E_v = par.Ev / Constants::Hartree - E1;
    else
      result.E_v = par.Ev / Constants::Hartree - E2;
  }



  //conduction-valence band coupling. may be needed only for optics
  result.P1 = std::sqrt(par.Ep_1 * 0.5 / Constants::Hartree);
  result.P2 = std::sqrt(par.Ep_2 * 0.5 / Constants::Hartree);


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
