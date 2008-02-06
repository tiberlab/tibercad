#include "ZbSemiconductor.h"
#include "getpot.h"
#include "Alloy.h"
#include "Database.h"
 
using namespace std;

typedef std::complex<double> Complex;
extern "C" 
{ 
  //ZHEEV( JOBZ, UPLO,  N, A,           LDA, W, WORK, LWORK, RWORK,INFO )
  void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& N, double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info);
};

const double ZbSemiconductor::Hartree = 27.2113961;
//---------------------------------------------//

ZbSemiconductor::ZbDDparameters& ZbSemiconductor::get_parameters()
{
  return(par);
}


//--------------------------------------------//
void ZbSemiconductor::do_init()
{
  Semiconductor::do_init();

  ModelOptions & options = get_options ();
  
  //-----------------------------------------------------------------------------
  //parameters
  {
    par.EgGamma = get_parameter("Eg_G",  par.EgGamma );
    par.EgL     = get_parameter("Eg_L",  par.EgL);
    par.EgX     = get_parameter("Eg_X",  par.EgX);
    par.Ev      = get_parameter("E_v",   par.Ev);
  
    par.m_G     = get_parameter("m_G",    par.m_G);
    par.m_l_L   = get_parameter("m_L_l",  par.m_l_L);
    par.m_t_L   = get_parameter("m_L_t",  par.m_t_L);
    par.m_l_X   = get_parameter("m_X_l",  par.m_l_X);
    par.m_t_X   = get_parameter("m_X_t",  par.m_t_X);
    
    par.a_c     = get_parameter("a_c", par.a_c);
    par.a_v     = get_parameter("a_v", par.a_v );
    par.b       = get_parameter("b",   par.b);
    par.d       = get_parameter("d",   par.d);
    
    par.delta  = get_parameter("delta",   par.delta);
    par.gamma1 = get_parameter("gamma1", par.gamma1);
    par.gamma2 = get_parameter("gamma2", par.gamma2);
    par.gamma3 = get_parameter("gamma3", par.gamma3);
    
    par.def_vol_X   = get_parameter("abs_def_pot_X",      par.def_vol_X );
    par.def_uniax_X = get_parameter("uniax_def_pot_X",  par.def_uniax_X);
    par.def_vol_L   = get_parameter("abs_def_pot_L",      par.def_vol_L);
    par.def_uniax_L = get_parameter("uniax_def_pot_L",  par.def_uniax_L);

    par.Ep = options.get_option("Ep", par.Ep);
  }
  //------------------------------------------------------------------------------
  //bowing
  {
    bow.EgGamma = get_parameter("bow_Eg_G",  bow.EgGamma );
    bow.EgL     = options.get_option("bow_Eg_L",  bow.EgL);
    bow.EgX     = options.get_option("bow_Eg_X",  bow.EgX);
    bow.Ev      = get_parameter("bow_E_v",   bow.Ev);
    
    bow.m_G     = options.get_option("bow_m_G",    bow.m_G);
    bow.m_l_L   = options.get_option("bow_m_L_l",  bow.m_l_L);
    bow.m_t_L   = options.get_option("bow_m_L_t",  bow.m_t_L);
    bow.m_l_X   = options.get_option("bow_m_X_l",  bow.m_l_X);
    bow.m_t_X   = options.get_option("bow_m_X_t",  bow.m_t_X);
    
    bow.a_c     = options.get_option("bow_a_c", bow.a_c);
    bow.a_v     = options.get_option("bow_a_v", bow.a_v );
    bow.b       = options.get_option("bow_b",   bow.b);
    bow.d       = options.get_option("bow_d",   bow.d);
    
    bow.delta  = options.get_option("bow_delta",   bow.delta);
    bow.gamma1 = options.get_option("bow_gamma1", bow.gamma1);
    bow.gamma2 = options.get_option("bow_gamma2", bow.gamma2);
    bow.gamma3 = options.get_option("bow_gamma3", bow.gamma3);
    
    bow.def_vol_X   = options.get_option("abs_def_pot_X",      bow.def_vol_X );
    bow.def_uniax_X = options.get_option("uniax_def_pot_X",  bow.def_uniax_X);
    bow.def_vol_L   = options.get_option("abs_def_pot_L",      bow.def_vol_L);
    bow.def_uniax_L = options.get_option("uniax_def_pot_L",  bow.def_uniax_L);

    bow.Ep = options.get_option("bow_Ep_1", bow.Ep);
  }


  
}



//-------------------------------------------//
void ZbSemiconductor::read_database( )
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());


  // defaults for GaAs
  par.EgGamma = data("Eg_G", 1.519);
  par.EgL = data("Eg_L", 1.815);
  par.EgX = data("Eg_X", 1.981);
  par.Ev = data("E_v", 1.346);
  
  par.m_G = data("m_G", 0.067);
  par.m_l_L = data("m_L_l", 1.9);
  par.m_t_L = data("m_L_t", 0.0754);
  par.m_l_X = data("m_X_l", 1.3);
  par.m_t_X = data("m_X_t", 0.23);
  
  par.a_c = data("a_c", -9.36);
  par.a_v = data("a_v", -1.21);
  par.b = data("b", -2.0);
  par.d = data("d", -4.8);
  
  par.delta = data("delta", 0.341);
  par.gamma1 = data("gamma1", 6.98);
  par.gamma2 = data("gamma2", 2.06);
  par.gamma3 = data("gamma3", 2.93);
  
  par.def_vol_X = data("abs_def_pot_X", -0.16);
  par.def_uniax_X = data("uniax_def_pot_X", 14.26);
  par.def_vol_L = data("abs_def_pot_L", -4.91);
  par.def_uniax_L = data("uniax_def_pot_L", 6.5);

  par.Ep = data("Ep", 25.0);
  
}
//----------------------------------------------//
void ZbSemiconductor::copy_from (const PhysicalModelInterface *rhs)
{
 

  const ZbSemiconductor* mod = dynamic_cast<const ZbSemiconductor*> (rhs);
  
  par = mod->par;

  bow = mod->bow;

}


//--------------------------------------------//
void ZbSemiconductor::read_bowing_parameters()
{
 
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());


 
  bow.EgGamma = data("bow_Eg_G",0.0);
  bow.EgL = data("bow_Eg_L", 0.0);
  bow.EgX = data("bow_Eg_X", 0.0);
  bow.Ev = data("bow_E_v", 0.0);
  
  bow.m_G = data("bow_m_G", 0.0);
  bow.m_l_L = data("bow_m_L_l", 0.0);
  bow.m_t_L = data("bow_m_L_t", 0.0);
  bow.m_l_X = data("bow_m_X_l", 0.0);
  bow.m_t_X = data("bow_m_X_t", 0.0);
  
  bow.a_c = data("bow_a_c", 0.0);
  bow.a_v = data("bow_a_v", 0.0);
  bow.b = data("bow_b", 0.0);
  bow.d = data("bow_d", 0.0);
  
  bow.delta = data("bow_delta", 0.0);
  bow.gamma1 = data("bow_gamma1", 0.0);
  bow.gamma2 = data("bow_gamma2", 0.0);
  bow.gamma3 = data("bow_gamma3", 0.0);
  
  bow.def_vol_X = data("bow_abs_def_pot_X", 0.0);
  bow.def_uniax_X = data("bow_uniax_def_pot_X", 0.0);
  bow.def_vol_L = data("bow_abs_def_pot_L", 0.0);
  bow.def_uniax_L = data("bow_uniax_def_pot_L", 0.0);

  bow.Ep = data("bow_Ep", 0.0);
}


//---------------------------------------------//
void ZbSemiconductor::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  const ZbSemiconductor* modA = dynamic_cast<const ZbSemiconductor*> (comp_A);
  const ZbSemiconductor* modB = dynamic_cast<const ZbSemiconductor*> (comp_B);

  par.EgGamma = alloy( modA->par.EgGamma, modB->par.EgGamma, xa, bow.EgGamma);
  par.EgL     = alloy( modA->par.EgL, modB->par.EgL, xa, bow.EgL);
  par.EgX     = alloy( modA->par.EgX, modB->par.EgX, xa, bow.EgX);
  par.Ev      = alloy( modA->par.Ev, modB->par.Ev, xa, bow.Ev);
  
  par.m_G     = alloy( modA->par.m_G, modB->par.m_G, xa, bow.m_G);
  par.m_l_L   = alloy( modA->par.m_l_L, modB->par.m_l_L, xa, bow.m_l_L);
  par.m_t_L   = alloy( modA->par.m_t_L, modB->par.m_t_L, xa, bow.m_t_L);
  par.m_l_X   = alloy( modA->par.m_l_X, modB->par.m_l_X, xa, bow.m_l_X);
  par.m_t_X   = alloy( modA->par.m_t_X, modB->par.m_t_X, xa, bow.m_t_X);
  
  par.a_c     = alloy( modA->par.a_c, modB->par.a_c, xa, bow.a_c);
  par.a_v     = alloy( modA->par.a_v, modB->par.a_v, xa, bow.a_v);
  par.b       = alloy( modA->par.b, modB->par.b, xa, bow.b);
  par.d       = alloy( modA->par.d, modB->par.d, xa, bow.d);
  
  par.delta  =  alloy( modA->par.delta, modB->par.delta, xa, bow.delta);
  par.gamma1 =  alloy( modA->par.gamma1, modB->par.gamma1, xa, bow.gamma1);
  par.gamma2 =  alloy( modA->par.gamma2, modB->par.gamma2, xa, bow.gamma2);
  par.gamma3 =  alloy( modA->par.gamma3, modB->par.gamma3, xa, bow.gamma3);
  
  par.def_vol_X   = alloy( modA->par.def_vol_X, modB->par.def_vol_X, xa, bow.def_vol_X);
  par.def_uniax_X = alloy( modA->par.def_uniax_X, modB->par.def_uniax_X , xa, bow.def_uniax_X);
  par.def_vol_L   = alloy( modA->par.def_vol_L, modB->par.def_vol_L, xa, bow.def_vol_L);
  par.def_uniax_L = alloy( modA->par.def_uniax_L, modB->par.def_uniax_L , xa, bow.def_uniax_L);
  


}

//--------------------------------------------//
ZbSemiconductor::ZbSemiconductor( )
{

 
}


//-----------------------------------------------------------//

KPparams ZbSemiconductor::calculate_8x8_kp_params (void )
{


  //we start from 6x6 kp parameters  
  KPparams  result = calculate_6x6_kp_params();

 
 
  //------------------------------------------------------------------
  //CONDUCTION BAND
  //we renormalize conduction band quadratic part to free electron mass
  result.s1 = 1.0;
  result.s2 = 1.0;

  double Ev_top = par.Ev + ((1.0/3.0) * par.delta );

  result.E_c =  (Ev_top + par.EgGamma)/Hartree;

  

  //--------------------------------------------------------------------

  double Ep; //Ep = P^2 * 2.0;
  
  double s = 1.0;

  Ep = s*(1.0/par.m_G - 1.0)*
    par.EgGamma * ( (par.EgGamma +  par.delta )/(par.EgGamma + 2.0/3.0 * par.delta ) ) / Hartree;



  result.P1 = std::sqrt( 0.5*Ep);
  result.P2 = std::sqrt( 0.5*Ep);
 

  //result.P1 =  result.P2 = 0;

  //rescale L and N

  double t =   0.5*Ep/( (par.EgGamma +  par.delta/3.0)/Hartree );

  

  result.L1 += t;
  result.L2 += t;

  result.N1 += t;
  result.N2 += t;
  

  result.N1_yx = result.M1;  result.N2_yx = result.N1_yx;
  result.N1_xy = result.N1 - result.N1_yx; result.N2_xy = result.N1_xy;

  

  return(result);
}

//-----------------------------------------------------------//
 
KPparams ZbSemiconductor::calculate_6x6_kp_params (void )

{
 
  KPparams  result;



  //-----------------------------------------------
 // Valence band k.p parameters:

    /*  

          L = \frac{1}{2} (-\gamma_1 - 4 \gamma_2 - 1)  ; 

          M = \frac{1}{2} ( 2\gamma_2 - \gamma_1  - 1 ) ;  

          N = -3\gamma_3;

	  N_{yx} = M; 

          N_{xy} = N -  N_{yx};
          
     */
        

  result.L1 = 0.5 * (- par.gamma1 - 4.0 * par.gamma2 - 1.0); result.L2 = result.L1;
 
  result.M1 = 0.5 * (2 * par.gamma2 - par.gamma1 - 1.0); result.M2 = result.M1;        
  result.M3 = result.M1;

  result.N1 = -3.0 * par.gamma3; result.N2 = -3.0 * par.gamma3;

  result.N1_yx = result.M1;  result.N2_yx = result.N1_yx;

  
  result.N1_xy = result.N1 - result.N1_yx; result.N2_xy = result.N1_xy;
  //-------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------//
  //  Valence band deformation potential:

      /*
        l  =  a_v + 2b; 
        m  =  a_v  - b; 
        n  =  \sqrt{3} d. 
      */

  result.l1s = (par.a_v + 2.0 * par.b)/Hartree; result.l2s = result.l1s;

  result.m1s = (par.a_v - par.b)/Hartree;  result.m2s = result.m1s; result.m3s = result.m1s;

  result.n1s = (sqrt(3.0) * par.d)/Hartree; result.n2s = result.n1s;  
  //------------------------------------------------------------------------------//
  //Conduction band deformation potential
  
  result.axs = par.a_c / Hartree;
  result.azs = result.axs;
  //------------------------------------------------------------------------------//
  //  Averaged valence band energy 

  
  

  result.E_v = par.Ev / Hartree;

  //------------------------------------------------------------------------------//

  //spin-orbit energy
  result.d1 = 0.0 ; //no crystal field splitting
  result.d2 = (par.delta/3) / Hartree;
  result.d3 = (par.delta/3) / Hartree;

  result.P1  = std::sqrt(par.Ep * 0.5 / Hartree);
  result.P2  = result.P1;

  return(result);
 
}

//=================================================================================//




//============================================================================================//
