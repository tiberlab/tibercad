#include "KPbulkHamiltonian.h"
#include "getpot.h"
#include "Alloy.h"
#include "PhysicalModelInterface.h"




using namespace std;
//==================================================================

void KPbulkHamiltonian::nullify_parameters(void)
{

  par.L1 = 0.0;
  par.L2 = 0.0;
  par.M1 = 0.0;
  par.M2 = 0.0;
  par.M3 = 0.0;
  par.N1 = 0.0;
  par.N2 = 0.0;    
  par.P1 = 0.0;
  par.P2 = 0.0;
  par.s1 = 0.0;
  par.s2 = 0.0;
  par.E_c = 0.0;
  par.E_v = 0.0;
  par.d1 = 0.0;
  par.d2 = 0.0;
  par.d3 = 0.0;
  par.N1_xy = 0.0;
  par.N1_yx = 0.0; 
  par.N2_xy = 0.0; 
  par.N2_yx = 0.0;
  par.l1s = 0.0;
  par.l2s = 0.0;
  par.n1s = 0.0;  
  par.n2s = 0.0;
  par.m1s = 0.0;
  par.m2s = 0.0;
  par.m3s = 0.0;
  par.axs = 0.0;
  par.azs = 0.0;

  kpVVtermSymmetric = false;
  kpCVtermSymmetric = true;

  strainM = Tensor2Sym(0);
}
//================================================================
KPbulkHamiltonian::~KPbulkHamiltonian()
{
  PhysicalModelInterface::destroy(semiconductor);
}


//================================================================

KPbulkHamiltonian::KPbulkHamiltonian( )
  : model_name("6x6"),
    semiconductor(NULL),
    band_min(2),
    band_max(7)
{
} 



//=================================================================

void KPbulkHamiltonian::do_init() 
{

  EFAbulkHamiltonian::do_init();
  
  const ModelOptions& opt =  get_options ();

  PhysicalModelInterface::destroy(semiconductor);


 
  semiconductor = Semiconductor::create( get_material() -> get_structure(), opt  );

  semiconductor->set_material(get_material());

  semiconductor->init();
  

  model_name = opt.get_option("kp_model","6x6");

  kpVVtermSymmetric = opt.get_option("kpVVtermSymmetric", false);

  kpCVtermSymmetric = opt.get_option("kpCVtermSymmetric", true);
  
 
  nullify_parameters();

  if (model_name == "6x6")
    {
      band_min = 2;
      band_max = 7;
     
    }
  else
    { 
      if (model_name == "8x8")
	{
	  band_min = 0;
	  band_max = 7;

	}
      else
	{
	  cerr << "Wrong kp model:  " << model_name << "\n";
	  exit(1);
	}
    }



 


  short i1 = 0;

  for (short i = band_min; i <= band_max; i++) 
    {
      kp_bands.push_back(i);
      kp_bands_map.insert(make_pair (i,i1) );
      i1++;
    }

  //nullify strain
  strainM = Tensor2Sym(0);

  //prepare k.p parameter
  par = semiconductor->calculate_kp_params (model_name);

  //calculate general Hamiltonian 
  calculate_Hamiltonian_gen();

  //apply k|| vector (even if it is zero-vector !!!)
  calculate_Hamiltonian_k_par();

  //-----------------------------------------------
  //calculate optical operator
  calculate_optical_operator();
  calculate_optical_operator_k_par();
  //-----------------------------------------------


}



//==================================================================//
void KPbulkHamiltonian::do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
  const KPbulkHamiltonian* modA = dynamic_cast<const KPbulkHamiltonian*> (comp_A);
  const KPbulkHamiltonian* modB = dynamic_cast<const KPbulkHamiltonian*> (comp_B);


  PhysicalModelInterface::destroy(semiconductor);
  semiconductor = static_cast<Semiconductor*>(modA->semiconductor->copy());
  assert(semiconductor != NULL);
  semiconductor->set_material(get_material());
  semiconductor->init_alloy(modA->semiconductor, modB->semiconductor, xa);

  strainM = modA->strainM;
  model_name = modA->model_name;
  kpVVtermSymmetric = modA->kpVVtermSymmetric;
  kpCVtermSymmetric = modA->kpCVtermSymmetric;
  band_min = modA->band_min;
  band_max = modA->band_max;
  set_k_vector(modA->k_vector);

  set_rotation_matrix();

  //prepare k.p parameter
  par = semiconductor->calculate_kp_params(model_name);

  //calculate general Hamiltonian 
  calculate_Hamiltonian_gen();

  //apply k|| vector (even if it is zero-vector !!!)
  calculate_Hamiltonian_k_par();
  

  //-----------------------------------------------
  //calculate optical operator

  calculate_optical_operator();

  calculate_optical_operator_k_par();
  //-----------------------------------------------
}

//==================================================================//
 
void KPbulkHamiltonian::calculate_Hamiltonian_gen(void)
{
  //=================================================================

  
 
  //Initialization
  //--------------------------------------------------
  //create matrix 8x8
  const vector<KPbulkHamiltonian::MatrixElement>  temp(8);

 
  Ham.resize(8,temp);
  //--------------------------------------------------

  for (short i = 0; i < 8; i++)
    for (short j = 0; j < 8; j++)
      {
	//--------constant
	Ham[i][j].constant =  Complex(0.0,0.0);

	//--------linear
	for (short i1 = 0; i1 < 3; i1++)  
	  {
	    Ham[i][j].linear_left[i1] =  Complex(0.0, 0.0);
	    Ham[i][j].linear_right[i1] =  Complex(0.0, 0.0);
	  }
        //--------quadratic
	for (short i1 = 0; i1 < 3; i1++)
	  for (short j1 = 0; j1 < 3; j1++)
	    Ham[i][j].quad[i1][j1] =  Complex(0.0, 0.0);

      }

 
  //=================================================================
  //-------H_vv + h^2/2m(kx^2+ky^2+kz^2) in crystal system-----------//
  Ham[2][2].quad[0][0]  = par.L1  + 0.5;
  Ham[2][2].quad[1][1]  = par.M1  + 0.5;
  Ham[2][2].quad[2][2]  = par.M2  + 0.5;

  Ham[3][3].quad[0][0] = par.M1 + 0.5;
  Ham[3][3].quad[1][1] = par.L1 + 0.5;
  Ham[3][3].quad[2][2] = par.M2 + 0.5;

  Ham[4][4].quad[0][0] = par.M3 + 0.5;
  Ham[4][4].quad[1][1] = par.M3 + 0.5;
  Ham[4][4].quad[2][2] = par.L2 + 0.5;
  
  if (kpVVtermSymmetric)
    {
      Ham[2][3].quad[0][1] = par.N1/2.0;  Ham[2][3].quad[1][0]= par.N1/2.0;
      Ham[2][4].quad[0][2] = par.N2/2.0;  Ham[2][4].quad[2][0]= par.N2/2.0;
      Ham[3][4].quad[1][2] = par.N2/2.0;  Ham[3][4].quad[2][1]= par.N2/2.0;
    }
  else
    {
      Ham[2][3].quad[0][1] = par.N1_xy;  Ham[2][3].quad[1][0]= par.N1_yx;
      Ham[2][4].quad[0][2] = par.N2_xy;  Ham[2][4].quad[2][0]= par.N2_yx;
      Ham[3][4].quad[1][2] = par.N2_xy;  Ham[3][4].quad[2][1]= par.N2_yx;
     
    }

  /*
   was in fortran
    Ham_quad(4,3,:,:) = TRANSPOSE(CONJG(Ham_quad(3,4,:,:)))
    Ham_quad(5,4,:,:) = TRANSPOSE(CONJG(Ham_quad(4,5,:,:)))
    Ham_quad(5,3,:,:) = TRANSPOSE(CONJG(Ham_quad(3,5,:,:)))
  */
  
  for (short i = 0; i< 3; i++)
     for (short j = 0; j< 3; j++)
       {
	 Ham[3][2].quad[i][j] = conj(Ham[2][3].quad[j][i]);
	 Ham[4][3].quad[i][j] = conj(Ham[3][4].quad[j][i]);
	 Ham[4][2].quad[i][j] = conj(Ham[2][4].quad[j][i]);
       }
  //=====================================================================!
  //--------Valence band strain ------------------------------------

 

  Ham[2][2].constant = par.l1s*strainM(1,1)+par.m1s*strainM(2,2)+par.m2s*strainM(3,3);
  Ham[3][3].constant = par.m1s*strainM(1,1)+par.l1s*strainM(2,2)+par.m2s*strainM(3,3); 
  Ham[4][4].constant = par.m3s*strainM(1,1)+par.m3s*strainM(2,2)+par.l2s*strainM(3,3); 
  Ham[2][3].constant = par.n1s*strainM(2,1) ; 
  Ham[3][2].constant = par.n1s*strainM(2,1) ;
  Ham[2][4].constant = par.n2s*strainM(3,1) ;  
  Ham[4][2].constant = par.n2s*strainM(3,1) ;
  Ham[3][4].constant = par.n2s*strainM(3,2) ;  
  Ham[4][3].constant = par.n2s*strainM(3,2) ;


 

  //=================================================================
  // Crystal field splitting
  Ham[2][2].constant +=  par.d1 ; 
  Ham[3][3].constant +=  par.d1 ; 
  //=================================================================
  
  //  !---------Valence part without spin-orbit-------------------------!

  /*
    Ham_quad(6:8,6:8,:,:)=Ham_quad(3:5,3:5,:,:)
    Ham_const(6:8,6:8)=Ham_const(3:5,3:5)
  */
  
  for (short i = 5 ; i <= 7; i++)
    for (short j = 5 ; j <= 7; j++)
      {
       
	Ham[i][j].constant = Ham[i-3][j-3].constant;

	for (short i1 = 0; i1 < 3; i1++)
	  for (short j1 = 0; j1 < 3; j1++)
	    Ham[i][j].quad[i1][j1] = Ham[i-3][j-3].quad[i1][j1];
      }

  //=======================================================================//
  //-------------H_cc------------------------------------------------------//
  Ham[0][0].quad[0][0]= 0.5*par.s2;       
  Ham[0][0].quad[1][1]= 0.5*par.s2;     
  Ham[0][0].quad[2][2]= 0.5*par.s1;

  Ham[1][1].quad[0][0]= 0.5*par.s2;       
  Ham[1][1].quad[1][1]= 0.5*par.s2;     
  Ham[1][1].quad[2][2]= 0.5*par.s1;


 
  //-----------------------------------------------------------------------//

 
  //cv part
  if (kpCVtermSymmetric)
    {
   

      Ham[0][2].linear_left[0] = par.P2 * 0.5 *  Complex(0.0, 1.0);
      Ham[0][3].linear_left[1] = par.P2 * 0.5 *  Complex(0.0, 1.0);
      Ham[0][4].linear_left[2] = par.P1 * 0.5 *  Complex(0.0, 1.0);

      Ham[0][2].linear_right[0] = par.P2* 0.5 *  Complex(0.0, 1.0);
      Ham[0][3].linear_right[1] = par.P2* 0.5 *  Complex(0.0, 1.0);
      Ham[0][4].linear_right[2] = par.P1* 0.5 *  Complex(0.0, 1.0);



      Ham[1][5].linear_left[0] = par.P2* 0.5 *  Complex(0.0, 1.0);
      Ham[1][6].linear_left[1] = par.P2* 0.5 *  Complex(0.0, 1.0);
      Ham[1][7].linear_left[2] = par.P1* 0.5 *  Complex(0.0, 1.0);

      Ham[1][5].linear_right[0] = par.P2* 0.5 *  Complex(0.0, 1.0);
      Ham[1][6].linear_right[1] = par.P2* 0.5 *  Complex(0.0, 1.0);
      Ham[1][7].linear_right[2] = par.P1* 0.5 *  Complex(0.0, 1.0);



     

      for (short i = 0; i < 2; i++)
	for (short j = 2; j < 8; j++)
	  for (short  i1 = 0; i1<=2; i1++)
	    {
	        Ham[j][i].linear_left[i1] =  conj(Ham[i][j].linear_left[i1]);
	        Ham[j][i].linear_right[i1] = conj(Ham[i][j].linear_right[i1]);

	     
	      
	    }



      
     
      


    }
  else
    {
      Ham[0][2].linear_left[0]= par.P2 *  Complex(0, 1.0);  Ham[1][5].linear_left[0]= par.P2 *  Complex(0, 1.0); 
      Ham[0][3].linear_left[1]= par.P2 *  Complex(0, 1.0);  Ham[1][6].linear_left[1]= par.P2 *  Complex(0, 1.0);
      Ham[0][4].linear_left[2]= par.P1 *  Complex(0, 1.0);  Ham[1][7].linear_left[1]= par.P1 *  Complex(0, 1.0);
     
      for (short i = 0; i< 2; i++)
	for (short j = 2; j< 8; j++)
	  for (short  i1 = 0; i1<=2; i1++)
	    Ham[j][i].linear_right[i1] = conj(Ham[i][j].linear_left[i1]);

    }
  //-----------------------------------------------------------------------//
  // absolute band edges 
  //

  //conduction band with strain
  Ham[0][0].constant +=  par.E_c;   Ham[1][1].constant +=  par.E_c;

  //valence band without strain
 

  for (short i = 2; i < 8; i++)
    { 
     
      Ham[i][i].constant += par.E_v; //E_v is averaged
     
    }
  //-----------------------------------------------------------------------//
  
  //------------------spin-orbit interaction------------------------//

  /*Ham_const(3,4) = Ham_const(3,4)  + (0d0,-1d0)*d2*/  Ham[2][3].constant   +=  Complex(0.0,-1.0)*par.d2; 
  /*Ham_const(4,3) = Ham_const(4,3)  + (0d0, 1d0)*d2*/  Ham[3][2].constant   +=  Complex(0.0, 1.0)*par.d2;

  /*Ham_const(3,8) = Ham_const(3,8)  +  d3*/     Ham[2][7].constant += par.d3;
  /*Ham_const(8,3) = Ham_const(8,3)  +  d3*/     Ham[7][2].constant += par.d3;

  /*Ham_const(6,5) = Ham_const(6,5)  -  d3*/     Ham[5][4].constant += -par.d3;
  /*Ham_const(5,6) = Ham_const(5,6)  -  d3*/     Ham[4][5].constant += -par.d3;

  /* Ham_const(7,5) = Ham_const(7,5)  + (0d0,-1d0)*d3*/ Ham[6][4].constant +=  Complex(0.0,-1.0)*par.d3;
  /* Ham_const(5,7) = Ham_const(5,7)  + (0d0, 1d0)*d3*/ Ham[4][6].constant +=  Complex(0.0,1.0)*par.d3;

  /* Ham_const(7,6) = Ham_const(7,6)  + (0d0,-1d0)*d2*/ Ham[6][5].constant +=  Complex(0.0,-1.0)*par.d2;
  /* Ham_const(6,7) = Ham_const(6,7)  + (0d0, 1d0)*d2*/ Ham[5][6].constant +=  Complex(0.0,1.0)*par.d2;

  /* Ham_const(8,4) = Ham_const(8,4)  + (0d0, 1d0)*d3*/ Ham[7][3].constant +=  Complex(0.0,1.0)*par.d3;
  /* Ham_const(4,8) = Ham_const(4,8)  + (0d0,-1d0)*d3*/ Ham[3][7].constant +=  Complex(0.0,-1.0)*par.d3;
  //-----------------------------------------------------------------!
  //Transformation from crystal to calculation system



  for (short i = 0; i < 8; i++)
  {  
    for (short j = 0; j < 8; j++)
    {
      rotate_linear( Ham[i][j].linear_left);
      rotate_linear( Ham[i][j].linear_right);
      rotate_quad(Ham[i][j].quad);
	
    }
  }

 
  //-----------------------------------------------------------------!

}

void KPbulkHamiltonian::calculate_optical_operator(void)
{
  //nullify P matrix
  P.resize(3);
  for (short i = 0 ; i < 3; i++)  
    {
      P[i].resize(8);
      for (short j = 0 ; j < 8; j++)	P[i][j].resize(8);
    }

  for (short pol = 0; pol < 3; pol++)
    for (short i = 0; i < 8; i++)
      for (short j = 0; j < 8; j++)
	{
	  //--------constant
	  P[pol][i][j].constant =  Complex(0.0,0.0);
	  
	  //--------linear
	  for (short i1 = 0; i1 < 3; i1++)  
	    {
	      P[pol][i][j].linear_left[i1] =  Complex(0.0, 0.0);
	      P[pol][i][j].linear_right[i1] =  Complex(0.0, 0.0);
	    }
	  
	  
	}

  
 


  //obtain P matrix
  //we calculate derivative of H(k) matrix
  for (short band1 = 0; band1 < 8; band1++)
    for (short band2 = 0; band2 < 8; band2++)
      {
	//derivative of the linear term 
	for (short polariz = 0; polariz < 3;  polariz++)
	  P[polariz][band1][band2].constant += Ham[band1][band2].linear_left[polariz] + Ham[band1][band2].linear_right[polariz];
	

	//derivative of the quadratic term
	for (short polariz = 0; polariz < 3;  polariz++)
	  for (short i1= 0; i1 < 3;  i1++)
	    {   
	      P[polariz][band1][band2].linear_left[i1]   += Ham[band1][band2].quad[polariz][i1];
	      P[polariz][band1][band2].linear_right[i1]  += Ham[band1][band2].quad[i1][polariz];
	    }

      }

  //----------------------------
  P_gen = P; //k|| = 0

}
  








//-------------------------------------------------------//

void KPbulkHamiltonian:: calculate_Hamiltonian_k_par (void)
{
  //allocation of the result
  //Initialization
   vector< vector<MatrixElement > > result = Ham;
  //--------------------------------------------------//


  for (short i = band_min; i <= band_max; i++)
    for (short j = band_min; j <= band_max; j++)
    {
	 
      //------we have to change constant term
      for (short i1 = 0; i1 < 3; i1++)
      {
	result[i][j].constant += Ham[i][j].linear_left[i1]  * k_vector[i1];
	result[i][j].constant += Ham[i][j].linear_right[i1] * k_vector[i1];
	for (short j1 = 0; j1 < 3; j1++)
	{
	  result[i][j].constant += Ham[i][j].quad[i1][j1] * k_vector[i1] * k_vector[j1]; 
	}
      }

      //------we have to change linear term
      
      for (short i1 = 0; i1 < 3; i1++)
	for (short j1 = 0; j1 < 3; j1++)
	{
	  result[i][j].linear_left[i1]  += Ham[i][j].quad[i1][j1] * k_vector[j1];
	  result[i][j].linear_right[j1] += Ham[i][j].quad[i1][j1] * k_vector[i1];
	}

    }
  

  //-------------------------------------------------//

  Hamiltonian.resize(band_max - band_min + 1);
  for (short i = 0; i <= band_max - band_min; i++)  Hamiltonian[i].resize(band_max - band_min + 1);


  for (short i = 0; i <= band_max - band_min; i++)
     for (short j = 0; j <= band_max - band_min; j++)
       Hamiltonian[i][j] = result[i + band_min][j + band_min];
  
  //-----------------------------------------------//
  Hamiltonian_without_strain_pot = Hamiltonian;
  

}

//-------------------------------------------------------//
void KPbulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{

  

  const vector<Complex>  strain_Ham_Bir_Pikus1(8,Complex(0.0,0.0));
  
  vector< vector<Complex > > strain_Ham_Bir_Pikus;
  strain_Ham_Bir_Pikus.resize(8,strain_Ham_Bir_Pikus1);


 
 
  

  //---------------------------------------------------------------
  //  strain 
  //---------------------------------------------------------------
  //valence band
  

  strain_Ham_Bir_Pikus[2][2]  = par.l1s*strain_crystal(1,1)+par.m1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3);
  strain_Ham_Bir_Pikus[3][3]  = par.m1s*strain_crystal(1,1)+par.l1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3); 
  strain_Ham_Bir_Pikus[4][4]  = par.m3s*strain_crystal(1,1)+par.m3s*strain_crystal(2,2)+par.l2s*strain_crystal(3,3); 

  strain_Ham_Bir_Pikus[2][3]  = par.n1s*strain_crystal(2,1) ; 
  strain_Ham_Bir_Pikus[3][2]  = par.n1s*strain_crystal(2,1) ;
  strain_Ham_Bir_Pikus[2][4]  = par.n2s*strain_crystal(3,1) ;  
  strain_Ham_Bir_Pikus[4][2]  = par.n2s*strain_crystal(3,1) ;
  strain_Ham_Bir_Pikus[3][4]  = par.n2s*strain_crystal(3,2) ;  
  strain_Ham_Bir_Pikus[4][3]  = par.n2s*strain_crystal(3,2) ;

  strain_Ham_Bir_Pikus[5][5]  = par.l1s*strain_crystal(1,1)+par.m1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3);
  strain_Ham_Bir_Pikus[6][6]  = par.m1s*strain_crystal(1,1)+par.l1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3); 
  strain_Ham_Bir_Pikus[7][7]  = par.m3s*strain_crystal(1,1)+par.m3s*strain_crystal(2,2)+par.l2s*strain_crystal(3,3); 

  strain_Ham_Bir_Pikus[5][6]  = par.n1s*strain_crystal(2,1) ; 
  strain_Ham_Bir_Pikus[6][5]  = par.n1s*strain_crystal(2,1) ;
  strain_Ham_Bir_Pikus[5][7]  = par.n2s*strain_crystal(3,1) ;  
  strain_Ham_Bir_Pikus[7][5]  = par.n2s*strain_crystal(3,1) ;
  strain_Ham_Bir_Pikus[6][7]  = par.n2s*strain_crystal(3,2) ;  
  strain_Ham_Bir_Pikus[7][6]  = par.n2s*strain_crystal(3,2) ;
  
  //conduction band
  strain_Ham_Bir_Pikus[0][0]  = par.axs * ( strain_crystal(1,1) + strain_crystal(2,2) ) + par.azs * strain_crystal(3,3);
  strain_Ham_Bir_Pikus[1][1]  = par.axs * ( strain_crystal(1,1) + strain_crystal(2,2) ) + par.azs * strain_crystal(3,3);
  

  //------------------------------------------------
  //potential
  //------------------------------------------------

  
  for (short i = 0; i < 8 ; i++)
  {
    strain_Ham_Bir_Pikus[i][i] -= el_potential/Hartree;
  }
  

  //--------------------------------------------------
  //correction of the final Hamiltonian
  //--------------------------------------------------

  for (short i = 0; i <= band_max - band_min; i++)
    for (short j = 0; j <= band_max - band_min; j++)
      Hamiltonian[i][j].constant = Hamiltonian_without_strain_pot[i][j].constant + strain_Ham_Bir_Pikus[i + band_min][j + band_min];
      

 

}
//======================================================//

//===============================================================/
void KPbulkHamiltonian::calculate_optical_operator_k_par(void)
{

  vector <vector< vector<MatrixElement > > > P = P_gen;

   for (short i = 0; i < 8; i++)
     for (short j = 0; j < 8; j++)
       {
	 //------we have to change constant term
	 for (short p = 0; p < 3; p++) 
	   for (short i1 = 0; i1 < 3; i1++)
	     {
	       P[p][i][j].constant += P_gen[p][i][j].linear_left[i1]  * k_vector[i1];
	       P[p][i][j].constant += P_gen[p][i][j].linear_right[i1] * k_vector[i1];
	       
	     }

       }

}


 



//-------------------------------------------------------//
const std::vector< std::vector <std::vector<EFAbulkHamiltonian::MatrixElement> > > & KPbulkHamiltonian::get_optical_operator() const
{
  return(P);
}
