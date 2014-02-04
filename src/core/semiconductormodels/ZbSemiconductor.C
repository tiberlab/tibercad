
// $Id$

#include "ZbSemiconductor.h"
#include "Database.h"
#include "Material.h"
#include "Constants.h"

#include "tensor_value.h"

using namespace std;


//--------------------------------------------//
void ZbSemiconductor::do_init()
{
  Semiconductor::do_init();

  //-----------------------------------------------------------------------------
  //parameters (do_init comes ofter read_database so input file can override)
  get_parameter("Eg_G",  par.EgGamma );
  get_parameter("Eg_L",  par.EgL);
  get_parameter("Eg_X",  par.EgX);
  get_parameter("E_v",   par.Ev);

  par.m_G     = get_option("m_G",    par.m_G);
  par.m_l_L   = get_option("m_L_l",  par.m_l_L);
  par.m_t_L   = get_option("m_L_t",  par.m_t_L);
  par.m_l_X   = get_option("m_X_l",  par.m_l_X);
  par.m_t_X   = get_option("m_X_t",  par.m_t_X);

  par.a_c     = get_option("a_c", par.a_c);
  par.a_v     = get_option("a_v", par.a_v );
  par.b       = get_option("b",   par.b);
  par.d       = get_option("d",   par.d);

  par.delta  = get_option("delta",   par.delta);
  par.gamma1 = get_option("gamma1", par.gamma1);
  par.gamma2 = get_option("gamma2", par.gamma2);
  par.gamma3 = get_option("gamma3", par.gamma3);

  if (_kp_model == "2x2")
  {
    double m_v = get_option("m_v", 0.0);
    if (m_v != 0.0)
    {
      par.gamma1 = 1.0 / m_v;
      par.gamma2 = 0.0;
    }

    double m_c = get_option("m_c", 0.0);
    if (m_c != 0.0)
    {
      par.m_G = m_c;
    }
  }

  par.def_vol_X   = get_option("abs_def_pot_X",      par.def_vol_X );
  par.def_uniax_X = get_option("uniax_def_pot_X",  par.def_uniax_X);
  par.def_vol_L   = get_option("abs_def_pot_L",      par.def_vol_L);
  par.def_uniax_L = get_option("uniax_def_pot_L",  par.def_uniax_L);

  par.Ep = get_option("Ep", par.Ep);



  par.varshni_alpha_G = get_option("varshni_alpha_G", par.varshni_alpha_G);
  par.varshni_alpha_L = get_option("varshni_alpha_L", par.varshni_alpha_L);
  par.varshni_alpha_X = get_option("varshni_alpha_X", par.varshni_alpha_X);


  par.varshni_beta_G = get_option("varshni_beta_G",  par.varshni_beta_G);
  par.varshni_beta_L = get_option("varshni_beta_L",  par.varshni_beta_L);
  par.varshni_beta_X = get_option("varshni_beta_X",  par.varshni_beta_X);


  if (!_couple_bands) 
  {
    _spurious = "none";
    par.Ep = 0.0;
  }

  //here zero temperature and work parameters coincide
  par_initial = par;
  

}



//-------------------------------------------//
void ZbSemiconductor::read_database( )
{

  const Database& db = get_database();

  // defaults for GaAs
  db.set_section("bandgap");
  par.EgGamma = db.get("Eg_G", 1.00);
  par.EgL = db.get("Eg_L", 1.00);
  par.EgX = db.get("Eg_X", 1.00);


  if (_consider_temperature)
  {
     par.varshni_alpha_G = db.get("varshni_alpha_G", 0.0);
     par.varshni_alpha_L = db.get("varshni_alpha_L", 0.0);
     par.varshni_alpha_X = db.get("varshni_alpha_X", 0.0);
     par.varshni_beta_G = db.get("varshni_beta_G", 0.0);
     par.varshni_beta_L = db.get("varshni_beta_L", 0.0);
     par.varshni_beta_X = db.get("varshni_beta_X", 0.0);
  }
  else
  {
    par.varshni_alpha_G = 0.0;
    par.varshni_alpha_L = 0.0;
    par.varshni_alpha_X = 0.0;
    par.varshni_beta_G = 0.0;
    par.varshni_beta_L = 0.0;
    par.varshni_beta_X = 0.0;
  }


  db.set_section("valenceband");
  par.Ev = db.get("E_v", 0.0);


  db.set_section("conductionband");
  
  RealTensor mten;
  db.get("m_G",mten,true);
  par.m_G = (mten(0,0)+mten(1,1)+mten(2,2))/3.0;

  par.m_l_L = db.get("m_L_l", 0.0);
  par.m_t_L = db.get("m_L_t", 0.0);
  par.m_l_X = db.get("m_X_l", 0.0);
  par.m_t_X = db.get("m_X_t", 0.0);


  // either use old 'kdotp' block or override with 'kdotp_6x6' 	  
  db.set_section("kdotp");
  par.delta = db.get("delta", 0.0);
  par.gamma1 = db.get("gamma1", 0.0);
  par.gamma2 = db.get("gamma2", 0.0);
  par.gamma3 = db.get("gamma3", 0.0);
  par.Ep = db.get("Ep", 0.0);
  par.Ep1 = db.get("Ep1", 0.00);
  par.Ep2 = db.get("Ep2", 0.00);
  par.m_G = db.get("m_c", par.m_G);
  par.m_c2 = db.get("m_c2", 1.0);
  par.Eg1 = db.get("Eg1", 0.00);

  //override with specific section e.g., 'kdotp_8x8'
  
  db.set_section("kdotp_"+_kp_model);

  par.delta = db.get("delta", par.delta);
  par.gamma1 = db.get("gamma1", par.gamma1);
  par.gamma2 = db.get("gamma2", par.gamma2);
  par.gamma3 = db.get("gamma3", par.gamma3);
  par.Ep = db.get("Ep", par.Ep);
  par.Ep1 = db.get("Ep1", par.Ep1);
  par.Ep2 = db.get("Ep2",par.Ep2);
  
  par.Eg1 = db.get("Eg1", par.Eg1);
  par.delta_c = db.get("delta_c", par.delta_c);
  par.delta_cf = db.get("delta_cf", par.delta_cf);

  par.m_G = db.get("m_c", par.m_G);
  par.m_c2 = db.get("m_c2", par.m_c2);
  

  db.set_section("deformation_potentials");
  par.a_c = db.get("a_c", 0.0);
  par.a_v = db.get("a_v", 0.0);
  par.b = db.get("b", 0.0);
  par.d = db.get("d", 0.0);

  par.def_vol_X = db.get("abs_def_pot_X", 0.0);
  par.def_uniax_X = db.get("uniax_def_pot_X", 0.0);
  par.def_vol_L = db.get("abs_def_pot_L", 0.0);
  par.def_uniax_L = db.get("uniax_def_pot_L", 0.0);

}


//--------------------------------------------//
void ZbSemiconductor::read_database_alloy()
{

  const Database& db = get_database();


  db.set_section("bandgap");
  bow.EgGamma = db.get("bow_Eg_G",0.0);
  bow.EgL = db.get("bow_Eg_L", 0.0);
  bow.EgX = db.get("bow_Eg_X", 0.0);

}




//--------------------------------------------//
ZbSemiconductor::ZbSemiconductor(const ModelOptions& options)
 : Semiconductor(options)
{


}

void ZbSemiconductor::do_calculate_kp_params (KPparams& par)
{
  if(_kp_model=="2x2") calculate_2x2_kp_params(par);
  else if(_kp_model=="6x6") calculate_6x6_kp_params(par);
  else if(_kp_model=="8x8") calculate_8x8_kp_params(par);
  else if(_kp_model=="14x14") calculate_14x14_kp_params(par);
}


void ZbSemiconductor::calculate_2x2_kp_params (KPparams&  result)
{
  calculate_6x6_kp_params(result);

  result.E_v = par.Ev / Constants::Hartree;

  result.s1 = 1.0/par.m_G;
  result.s2 = result.s1;

  result.E_c = (par.Ev + par.EgGamma)/Constants::Hartree;
  result.l1s = par.a_v/Constants::Hartree; result.l2s = result.l1s;
  result.m1s = par.a_v/Constants::Hartree; result.m2s = result.m1s;

  if (_spurious=="none")
  {
    result.s1 = 1.0/par.m_G - par.Ep/par.EgGamma;
    result.s2 = result.s1;
  }
  else if (_spurious=="Foreman")
  {
    result.s1 = 0.0;
    result.s2 = 0.0;
    par.Ep = par.EgGamma / par.m_G;
  }  
  else if (_spurious=="Chuang")
  {
    result.s1 = 1.0/par.m_c2;
    result.s2 = 1.0/par.m_c2;
    par.Ep = (1.0/par.m_G - result.s1)* par.EgGamma;
  }

  result.P1  = std::sqrt(par.Ep * 0.5 / Constants::Hartree);
  result.P2  = result.P1;
  
}
//=================================================================================//

void ZbSemiconductor::calculate_6x6_kp_params (KPparams&  result)
{

  //-----------------------------------------------
  // Valence band k.p parameters:

    /*
          L = \frac{1}{2} (-\gamma_1 - 4 \gamma_2 - 1)  ;

          M = \frac{1}{2} ( 2\gamma_2 - \gamma_1  - 1 ) ;

          N = -3\gamma_3;

	  N_{yx} = M;

          N_{xy} = N -  N_{yx};
     */

  result.L1 = 0.5 * (- par.gamma1 - 4.0 * par.gamma2 - 1.0);
  result.L2 = result.L1;
  result.M1 = 0.5 * (2 * par.gamma2 - par.gamma1 - 1.0);
  result.M2 = result.M1;
  result.M3 = result.M1;

  result.N1 = -3.0 * par.gamma3;
  result.N2 = -3.0 * par.gamma3;

  // Burd format operator ordering (see also Foreman)
  result.N1_yx = result.M1;
  result.N2_yx = result.N1_yx;
  result.N1_xy = result.N1 - result.N1_yx;
  result.N2_xy = result.N1_xy;


  //-------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------//
  //  Valence band deformation potential:
  //
  //    l  =  a_v + 2b;
  //    m  =  a_v  - b;
  //    n  =  \sqrt{3} d.
  //
  result.l1s = (par.a_v + 2.0 * par.b)/Constants::Hartree; result.l2s = result.l1s;

  result.m1s = (par.a_v - par.b)/Constants::Hartree;  result.m2s = result.m1s; result.m3s = result.m1s;

  result.n1s = (sqrt(3.0) * par.d)/Constants::Hartree; result.n2s = result.n1s;


  //------------------------------------------------------------------------------//
  //Conduction band deformation potential
  //
  result.axs = par.a_c / Constants::Hartree;
  result.azs = result.axs;


  //------------------------------------------------------------------------------//
  // valence bands center of mass
  //
  result.E_v = (par.Ev - par.delta / 3.0) / Constants::Hartree;


  //------------------------------------------------------------------------------//
  //spin-orbit energy
  result.d1 = 0.0 ; //no crystal field splitting
  result.d2 = (par.delta/3.0) / Constants::Hartree;
  result.d3 = (par.delta/3.0) / Constants::Hartree;

  result.P1  = std::sqrt(par.Ep * 0.5 / Constants::Hartree);
  result.P2  = result.P1;


}

//=================================================================================//
void ZbSemiconductor::calculate_8x8_kp_params (KPparams&  result )
{

  //we start from 6x6 kp parameters
  calculate_6x6_kp_params(result);

  //------------------------------------------------------------------
  // CONDUCTION BAND
  // Treatment of the spurious solutions
  //
  // No special care. Normal k.p model
  //
  if (_spurious=="none")
  {
    result.s1 = 1.0/par.m_G - par.Ep/3.0 * (2.0/par.EgGamma + 1.0/(par.EgGamma + par.delta) );
    
    result.s2 = result.s1;
  }
  else if (_spurious=="Foreman")
  {
    // Foreman: Ac = 0.0
    result.s1 = 0.0;
    result.s2 = 0.0;

    par.Ep = 3.0/par.m_G * (par.EgGamma +  par.delta )/(3.0 * par.EgGamma + 2.0 * par.delta ) 
              * par.EgGamma;
    
  }
  else if (_spurious=="Chuang")
  {
    // Chuang/Povolotskty: Ac = 1/2.0  (division by 2 is done on assembly)
    result.s1 = 1.0/par.m_c2;
    result.s2 = 1.0/par.m_c2;
       
    par.Ep = 3.0*(1.0/par.m_G - result.s1)*(par.EgGamma+par.delta)/(3.0*par.EgGamma + 2.0*par.delta )
              * par.EgGamma;
    
  }
  result.P1  = std::sqrt(par.Ep * 0.5 / Constants::Hartree);
  result.P2 = result.P1;
 
  //if (_spurious == "none") 
  //{
  //  std::cout<<"m_c= "<<1.0/result.s1<<std::endl;
  // }
  //else
  //{
  //  std::cout<<"Ep= "<<par.Ep<<std::endl;
  //}  
  //--------------------------------------------------------------------
  //rescale L and N
  // this comes from gamma1= gamma1_L - 1/3 Ep / (Eg + delta/3)
  //                 gamma2= gamma2_L - 1/6 Ep / (Eg + delta/3) = gamma3
  double t = 0.5*par.Ep/( (par.EgGamma +  par.delta/3.0));

  result.L1 += t;
  result.L2 += t;


  result.N1 += t;
  result.N2 += t;

  // Burd format operator ordering
  result.N1_yx = result.M1;                  //N-
  result.N2_yx = result.N1_yx;
  result.N1_xy = result.N1 - result.N1_yx;   //N+
  result.N2_xy = result.N1_xy;


  // from database we get the valence band maximum (VBM)

  result.E_c =  (par.Ev + par.EgGamma)/Constants::Hartree;

}

//=================================================================================//
void ZbSemiconductor::calculate_14x14_kp_params(KPparams&  result)
{

  //we start from 8x8 kp parameters
  calculate_8x8_kp_params(result);

  // 14x14 correction
   double r = par.Ep1/3.0 * ( 1.0/(par.Eg1-par.EgGamma) + 2.0/(par.Eg1-par.EgGamma+par.delta_c) );
  //double r = 0.0;

  if (_spurious=="none")
  {
      result.s1 += r;
      result.s2 = result.s1;
  }
  else if (_spurious=="Foreman")
  {
    // Foreman: Ac = 0.0
    result.s1 = 0.0;
    result.s2 = 0.0;

    par.Ep = 3.0*(1.0/par.m_G + r) * (par.EgGamma +  par.delta )/(3.0 * par.EgGamma + 2.0 * par.delta ) 
              * par.EgGamma;
    
  }
  else if (_spurious=="Chuang")
  {
    result.s1 = 1.0/par.m_c2;
    result.s2 = 1.0/par.m_c2;
       
    par.Ep = 3.0*(1.0/par.m_G - result.s1 + r)*(par.EgGamma+par.delta)/(3.0*par.EgGamma + 2.0*par.delta )
              * par.EgGamma;
    
  }

  if (_spurious == "none") 
  {
    std::cout<<"Ac= "<<result.s1<<std::endl;
  }
  else
  {
    std::cout<<"Ep= "<<par.Ep<<std::endl;
  }  
  //rescale L and N
  // this comes from gamma1_14x14= gamma1_8x8 - 2/3 Ep2 / (Eg1 + delta/3 + 2 delta' /3)
  //                 gamma2_14x14= gamma2_8x8 + 1/6 Ep2 / (Eg1 + delta/3 + 2 delta' /3)
  //                 gamma3_14x14= gamma3_8x8 - 1/6 Ep2 / (Eg1 + delta/3 + 2 delta' /3)
  //
  double q = 0.5* par.Ep2/( (par.Eg1 +  par.delta/3.0 + 2.0*par.delta_c/3.0));
  //double q =0.0;

  result.M1 += q;
  result.M2 = result.M1;
  result.M3 = result.M1;

  result.N1 += q;
  result.N2 += q;

  // Burd format operator ordering
  result.N1_yx = result.M1;
  result.N2_yx = result.N1_yx;
  result.N1_xy = result.N1 - result.N1_yx;
  result.N2_xy = result.N1_xy;

  // restore standard P0
  result.P1 = std::sqrt(0.5*par.Ep/Constants::Hartree);
  result.P2 = result.P1;

  result.P1_c = std::sqrt(0.5*par.Ep1/Constants::Hartree);
  result.P2_c = std::sqrt(0.5*par.Ep2/Constants::Hartree);

  result.E_c1 = (par.Ev + par.Eg1 + 2.0/3.0*par.delta_c )/Constants::Hartree;


  result.d4 = par.delta_c/3.0 / Constants::Hartree;
  result.d5 = par.delta_cf/3.0 / Constants::Hartree;

  result.s3 = 0.0; //1.0/par.m_c2;
  result.s4 = 0.0; //1.0/par.m_c2;
  

}


//=================================================================================//
void ZbSemiconductor::apply_temperature(void)
{
  par = par_initial;



  if (get_material()->is_alloy())
  {



    const ZbSemiconductor::ZbDDparameters& parA_zero =  (dynamic_cast<const ZbSemiconductor*> (modelA))->get_initial_parameters();

    const ZbSemiconductor::ZbDDparameters& parB_zero =  (dynamic_cast<const ZbSemiconductor*> (modelB))->get_initial_parameters();



    double EgGammaA = parA_zero.EgGamma
      - parA_zero.varshni_alpha_G * _temperature*_temperature/(_temperature + parA_zero.varshni_beta_G);
    double EgGammaB = parB_zero.EgGamma
      - parB_zero.varshni_alpha_G * _temperature*_temperature/(_temperature + parB_zero.varshni_beta_G);

    par.EgGamma = alloy( EgGammaA, EgGammaB, _xa, bow.EgGamma);


    double EgXA = parA_zero.EgX
      - parA_zero.varshni_alpha_G * _temperature*_temperature/(_temperature + parA_zero.varshni_beta_G);
    double EgXB = parB_zero.EgX
      - parB_zero.varshni_alpha_G * _temperature*_temperature/(_temperature + parB_zero.varshni_beta_G);

    par.EgX = alloy( EgXA, EgXB, _xa, bow.EgX);



    double EgLA = parA_zero.EgL
      - parA_zero.varshni_alpha_G * _temperature*_temperature/(_temperature + parA_zero.varshni_beta_G);
    double EgLB = parB_zero.EgL
      - parB_zero.varshni_alpha_G * _temperature*_temperature/(_temperature + parB_zero.varshni_beta_G);

    par.EgL = alloy( EgLA, EgLB, _xa, bow.EgL);



  }
  else
  {



    par.EgGamma = par_initial.EgGamma - par.varshni_alpha_G * _temperature*_temperature/(_temperature + par.varshni_beta_G) ;
    par.EgX = par_initial.EgX - par.varshni_alpha_X * _temperature*_temperature/(_temperature + par.varshni_beta_X) ;
    par.EgL = par_initial.EgL - par.varshni_alpha_L * _temperature*_temperature/(_temperature + par.varshni_beta_L) ;
  }
}

//============================================================================================//
