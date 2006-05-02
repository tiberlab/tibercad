using namespace std;
#include "ZbDDsemiconductor.h"
typedef std::complex<double> Complex;
extern "C" 
{ 
  //ZHEEV( JOBZ, UPLO,  N, A,           LDA, W, WORK, LWORK, RWORK,INFO )
  void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& N, double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info);
};
const double ZbDDsemiconductor::Hartree;
//---------------------------------------------//

ZbDDsemiconductor::ZbDDsemiconductor(void): DDsemiconductor()
{
  par.Ev = 0;
  par.EgGamma = 0;
  par.EgL = 0;
  par.EgX = 0;
  par.gamma1 = 0;
  par.gamma2 = 0;
  par.gamma3 = 0; 
  par.m_t_L = 0;
  par.m_l_L = 0;    
  par.m_t_X = 0;
  par.m_l_X = 0;  
  par.m_G = 0;             
  par.a_c = 0;
  par.a_v = 0;
  par.b = 0;
  par.d = 0; 
  par.def_vol_X = 0; 
  par.def_uniax_X = 0; //volume and uniax deformation potential for X point
  par.def_vol_L = 0;
  par.def_uniax_L = 0; //volume and uniax deformation potential for L point
  par.delta = 0;


  

} 

//--------------------------------------------//
ZbDDsemiconductor::ZbDDsemiconductor(const ZbDDparameters& params): DDsemiconductor()
{

 

  par = params;

 
 
}

//---------------------------------------------//

void ZbDDsemiconductor::set_Eg(double EgGamma_1, double EgL_1, double EgX_1)
{
  par.EgGamma = EgGamma_1;
  par.EgL = EgL_1;
  par.EgX = EgX_1;

}

//--------------------------------------------//

void ZbDDsemiconductor::set_mass_Gamma(double m)
{
  par.m_G = m;
}

//--------------------------------------------//

void ZbDDsemiconductor::set_masses_L(double m_t, double m_l)
{
  par.m_t_L = m_t;
  par.m_l_L = m_l;
}

//-------------------------------------------//
void ZbDDsemiconductor::set_masses_X(double m_t, double m_l)
{
  par.m_t_X = m_t;
  par.m_l_X = m_l;
}

//-------------------------------------------//

void ZbDDsemiconductor::set_6x6kp_params(double gamma1_, double gamma2_, double gamma3_, double delta_)
{
  par.gamma1 = gamma1_;
  par.gamma2 = gamma2_;
  par.gamma3 = gamma3_;
  par.delta  = delta_;
}

//--------------------------------------------//
void ZbDDsemiconductor::set_deformation_parameters(double a_c_, double a_v_, double b_, double d_)
{
  par.a_c = a_c_;
  par.a_v = a_v_;
  par.b = b_;
  par.d = d_;
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
      Ec_L  += par.def_vol_L * trace(strain);
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
      band_ex.degeneracy = 16   ; /* spin degeneracy and  8 equivalent minima (
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
      Ec_X  += par.def_vol_X * trace(strain);
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

//-----------------------------------------------------------//


void ZbDDsemiconductor::set_Ev(const double Ev)
{
  par.Ev = Ev;
}


//-----------------------------------------------------------//
 
KPbulkHamiltonian::KPparams ZbDDsemiconductor::calculate_6x6_kp_params (void )

{
 
  KPbulkHamiltonian::KPparams  result;



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
  //------------------------------------------------------------------------------//
  //  Averaged valence band energy 

  
  

  result.E_v = par.Ev / Hartree;

  //------------------------------------------------------------------------------//

  //spin-orbit energy
  result.d1 = 0.0 ; //no crystal field splitting
  result.d2 = (par.delta/3) / Hartree;
  result.d3 = (par.delta/3) / Hartree;

  return(result);
 
}

//=================================================================================//
/*
vector< vector<double> >  ZbDDsemiconductor::calculate_vb_bulk_states(vector<Tensor1>& k_vector)
{
 

  vector< vector<double> > result;

  KPbulkHamiltonian  bulk;
  bulk.strainM = strain;

  //cerr << setw(16) << bulk.strainM ;
 
  KPbulkHamiltonian::KPparams params_kp ;
 
  params_kp =  calculate_6x6_kp_params();
  
  
  bulk.kpVVtermSymmetric = true;


  bulk.set_parameters( params_kp );
  Tensor2Gen RotMatrix;
  RotMatrix =  Tensor2Gen(1);

  bulk.set_rotation_matrix( RotMatrix );

  bulk.calculate_8x8kp_ham();

  std::vector<std::vector<KPbulkHamiltonian::element> >    Ham1;

  double kvec[3];


  std::complex<double> ham6x6matrix[6*6];

  unsigned int N = k_vector.size();

  vector <double> eigvals_calculated(6);
  
  for (short i1 = 0; i1 < N; i1++ )

    {

      kvec[0] =  k_vector[i1](1);  kvec[1] =  k_vector[i1](2);  kvec[2] =  k_vector[i1](3);


      Ham1 = bulk.get_Hamiltonian_k_par(kvec) ;

 

      for (short i = 0; i < 6; i++)
	{
	  for (short j = 0; j < 6; j++)
	    {
	     
	      ham6x6matrix[i + j*6] = Ham1[i+2][j+2].constant;
	     
	    }
	}

    
      char jobs = 'N';
      char UPLO = 'U'; 
      int  N = 6;
      double eigvals[6];
      std::complex<double> WORK[11];
      int LWORK = 11;
      double RWORK[16];
      int info;
      
      


      zheev_(jobs, UPLO, N, ham6x6matrix, N, eigvals, WORK, LWORK, RWORK, info); 
      if (info !=0 ) exit(1);


      for (short i = 0; i < 6; i++)
	eigvals_calculated[i] = eigvals[i]*Hartree;

      result.push_back(eigvals_calculated);
      

    }


  return(result);
}
*/
//--------------------------------------------------------//

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
	  //cerr << "imass\n";
	  //cerr << setw(14) << imass << "\n";

	  imass.invariants(&temp1, &temp2,&imass_DOS);

	  //cerr << "imass_DOS  " << imass_DOS << "\n";
	  
	  extremum.mass_DOS = std::pow(1.0/imass_DOS,1.0/3.0);
	  result.push_back(extremum);
	}
    }
  
 
  valence_band = result;
}
//---------------------------------------------------------//


//============================================================================================//
