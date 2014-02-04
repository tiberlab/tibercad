// $Id$

#include "WzSemiconductor.h"
#include "Database.h"
#include "Material.h"
#include "Constants.h"
#include "tensor_value.h"

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

  get_parameter("Eg_G", par.EgGamma );
  get_parameter("E_v", par.Ev );
  
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
  par.varshni_beta_G  = get_option("varshni_beta_G",  par.varshni_beta_G);
  
   
  // we override with hole band masses for 2x2
  if(_kp_model=="2x2")
  {
    RealTensor mten(0.0);
    get_option("m_v",mten);

    if (mten(0,0) != 0.0)
      par.A2 = par.A4 = -0.5 / mten(0,0);

    if (mten(2,2) != 0.0)
      par.A1 = par.A3 = -0.5 / mten(2,2);

    mten(0,0) = 0.0; mten(2,2) = 0.0;

    get_option("m_c",mten);
    if (mten(0,0) != 0.0) par.m_c_xx = mten(0,0);
    if (mten(2,2) != 0.0) par.m_c_zz = mten(2,2);

  }

  if (!_couple_bands) 
  {
    _spurious = "none";
    par.Ep_1 = 0.0;
    par.Ep_2 = 0.0;
  }


  //here zero temperature and work parameters coincide
  par_initial = par;
}


//--------------------------------------------------//
void  WzSemiconductor::read_database_alloy(void)
{

  const Database& db = get_database();

  db.set_section("bandgap");
  bow.EgGamma = db.get("bow_Eg_G", 0.0);

}


//--------------------------------------------------

void WzSemiconductor::read_database( )
{
  const Database& db = get_database();

  db.set_section("bandgap");
  par.EgGamma = db.get("Eg_G", 3.51);

  par.varshni_alpha_G = db.get("varshni_alpha_G", 0.0);
  par.varshni_beta_G  = db.get("varshni_beta_G", 0.0);


  db.set_section("valenceband");
  par.Ev = db.get("E_v", 0.0);

  db.set_section("conductionband");

  //db.get("m_c_xx", 0.0);
  //db.get("m_c_zz", 0.0);  
  // Read as tensor:
  RealTensor mten;
  db.get("m_G",mten,true);
  par.m_c_xx = mten(0,0);
  par.m_c_zz = mten(2,2);


  // read default 'kdotp' block  	  
  db.set_section("kdotp");

  par.A1 = db.get("A1", 0.0, true);
  par.A2 = db.get("A2", 0.0, true);
  par.A3 = db.get("A3", 0.0, true);
  par.A4 = db.get("A4", 0.0, true);
  par.A5 = db.get("A5", 0.0, true);
  par.A6 = db.get("A6", 0.0, true);  
  par.delta_s = db.get("delta_s", 0.0);
  par.delta_cr = db.get("delta_cr", 0.0);  
  par.Ep_1 = db.get("Ep_1", 0.0);
  par.Ep_2 = db.get("Ep_2", 0.0);
  
  // read from specific section (e.g. kdotp_8x8)
  db.set_section("kdotp_"+_kp_model);

  par.A1 = db.get("A1", par.A1);
  par.A2 = db.get("A2", par.A2);
  par.A3 = db.get("A3", par.A3);
  par.A4 = db.get("A4", par.A4);
  par.A5 = db.get("A5", par.A5);
  par.A6 = db.get("A6", par.A6); 
  par.delta_s = db.get("delta_s", par.delta_s);
  par.delta_cr = db.get("delta_cr", par.delta_cr); 
  par.Ep_1 = db.get("Ep_1", par.Ep_1);
  par.Ep_2 = db.get("Ep_2", par.Ep_2);
  

  db.set_section("deformation_potentials");
  mten = 0.0;
  db.get("a_c",mten,true);
  par.a_x = mten(0,0);
  par.a_z = mten(2,2);
  //par.a_x = db.get("a_x",0.0);
  //par.a_z = db.get("a_z",0.0);

  par.D1 = db.get("D1", par.D1);
  par.D2 = db.get("D2", par.D2);
  par.D3 = db.get("D3", par.D3);
  par.D4 = db.get("D4", par.D4);
  par.D5 = db.get("D5", par.D5);
  par.D6 = db.get("D6", par.D6);
  
}

void WzSemiconductor::do_calculate_kp_params (KPparams& par)
{
  if(_kp_model=="2x2") calculate_2x2_kp_params(par);
  else if(_kp_model=="6x6") calculate_6x6_kp_params(par);
  else if(_kp_model=="8x8") calculate_8x8_kp_params(par);
  else if(_kp_model=="14x14") calculate_14x14_kp_params(par);
}



void WzSemiconductor::calculate_2x2_kp_params (KPparams&  result)
{
  // Valence band k.p parameters:

  result.L1 = 0.5 * (par.A4 + par.A2 - 1.0);

  //result.L2 = 0.5 * (par.A1 - 1.0);

  result.M1 = 0.5 * (par.A4 + par.A2 - 1.0);

  result.M2 = 0.5 * (par.A1 + par.A3 - 1.0);

  //result.M3 = 0.5 * (par.A2 - 1.0);

  /*
  result.N1 = par.A5;

  result.N2 = par.A6/sqrt(2.0);

  result.N1_yx = result.M1;  result.N2_yx = result.M2;

  result.N1_xy = result.N1 - result.N1_yx; result.N2_xy = result.N2 - result.N2_yx;
  */

  result.l1s = (par.D4 + par.D2)/Constants::Hartree;
  result.m1s = (par.D4 + par.D2)/Constants::Hartree;
  result.m2s = (par.D1 + par.D3)/Constants::Hartree;

  //result.m3s = par.D2/Constants::Hartree;

  //result.n1s = 2.0 * par.D5/Constants::Hartree; result.n2s = sqrt(2.0) * par.D6/Constants::Hartree;


  //------------------------------------------------------------------------------//
  // conduction band strain
  result.axs = par.a_x / Constants::Hartree;
  result.azs = par.a_z / Constants::Hartree;

  result.E_v = par.Ev / Constants::Hartree;


  // conduction band
  result.s1 = 1.0/par.m_c_zz;
  result.s2 = 1.0/par.m_c_xx;
  result.E_c =  (par.Ev + par.EgGamma)/Constants::Hartree;

  // CB-VB coupling

  if (_spurious == "none")
  {
    result.s1 -= par.Ep_1 / par.EgGamma;
    result.s2 -= par.Ep_2 / par.EgGamma;
  }
  else if (_spurious == "Foreman")
  {
    // Foreman: Ac = 0.0
    result.s1 = 0.0;
    result.s2 = 0.0;

    par.Ep_1 = par.m_c_zz * par.EgGamma;

    par.Ep_2 = par.m_c_xx * par.EgGamma;
  }
  else if (_spurious == "Chuang")
  {
    //we renormalize conduction band quadratic part to free electron mass
    // Chuang/Povolotskty: Ac = 1/2.0  (division by 2 is done on assembly)

    result.s1 = 1.0;
    result.s2 = 1.0;

    par.Ep_1 = (1.0/par.m_c_zz - result.s1) * par.EgGamma;

    par.Ep_2 = (1.0/par.m_c_xx - result.s2) * par.EgGamma;
  }

  result.P1 = std::sqrt(0.5 * par.Ep_1 / Constants::Hartree);
  result.P2 = std::sqrt(0.5 * par.Ep_2 / Constants::Hartree);

}


void WzSemiconductor::calculate_6x6_kp_params (KPparams&  result)
{

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
  result.N1_yx = result.M1;  result.N2_yx = result.M2;

  result.N1_xy = result.N1 - result.N1_yx; result.N2_xy = result.N2 - result.N2_yx;
  //-------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------//

  result.l1s = (par.D5 +  par.D4 + par.D2)/Constants::Hartree; result.l2s = par.D1/Constants::Hartree;

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
  result.d2 = (par.delta_s/3.0) / Constants::Hartree;
  result.d3 = (par.delta_s/3.0) / Constants::Hartree;


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



  //conduction-valence band coupling. maybe needed only for optics
  result.P1 = std::sqrt(par.Ep_1 * 0.5 / Constants::Hartree);
  result.P2 = std::sqrt(par.Ep_2 * 0.5 / Constants::Hartree);


}



//--------------------------------------------------------------------
void WzSemiconductor::calculate_8x8_kp_params (KPparams&  result )
{
  //we start from 6x6 kp parameters
  calculate_6x6_kp_params(result);

  //------------------------------------------------------------------
  // treatment of spurious solutions
  if (_spurious=="none")
  {
    result.s1 = 1.0/par.m_c_zz - par.Ep_1/3.0 * (2.0/par.EgGamma + 1.0/(par.EgGamma + par.delta_s) );
    
    result.s2 = 1.0/par.m_c_xx - par.Ep_2/3.0 * (2.0/par.EgGamma + 1.0/(par.EgGamma + par.delta_s) );
  }
  else if (_spurious=="Foreman")
  {
    // Foreman: Ac = 0.0
    result.s1 = 0.0;
    result.s2 = 0.0;

    par.Ep_1 = 3.0/par.m_c_zz * (par.EgGamma +  par.delta_s )/(3.0 * par.EgGamma + 2.0 * par.delta_s ) 
              * par.EgGamma;

    par.Ep_2 = 3.0/par.m_c_xx * (par.EgGamma +  par.delta_s )/(3.0 * par.EgGamma + 2.0 * par.delta_s ) 
              * par.EgGamma;
    
  }
  else if (_spurious=="Chuang")
  {
    //we renormalize conduction band quadratic part to free electron mass
    // Chuang/Povolotskty: Ac = 1/2.0  (division by 2 is done on assembly)
    
    result.s1 = 1.0;
    result.s2 = 1.0;
    
    par.Ep_1 = (1.0/par.m_c_zz - result.s1)*
      par.EgGamma * ( (par.EgGamma +  par.delta_s )/(par.EgGamma + 2.0/3.0 * par.delta_s ) );
    
    par.Ep_2 = (1.0/par.m_c_xx - result.s2)*
      par.EgGamma * ( (par.EgGamma +  par.delta_s )/(par.EgGamma + 2.0/3.0 * par.delta_s ) );
  }
  //-----------------------------------------------------------------

  result.P1 = std::sqrt(0.5 * par.Ep_1 / Constants::Hartree);
  result.P2 = std::sqrt(0.5 * par.Ep_2 / Constants::Hartree);


  //-----------------------------------------------------------------
  //to check !
  //rescale L and N

  double t1 =   0.5*par.Ep_1/( (par.EgGamma +  par.delta_s/3.0) );
  double t2 =   0.5*par.Ep_2/( (par.EgGamma +  par.delta_s/3.0) );

  result.L1 += t1;
  result.L2 += t2;

  result.N1 += t1;
  result.N2 += t2;


  result.N1_yx = result.M1;  result.N2_yx = result.N1_yx;
  result.N1_xy = result.N1 - result.N1_yx; result.N2_xy = result.N1_xy;

 
  // from database we get the valence band maximum (VBM)
  result.E_c =  (par.Ev + par.EgGamma)/Constants::Hartree;

}

void WzSemiconductor::calculate_14x14_kp_params (KPparams&  result )
{
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
