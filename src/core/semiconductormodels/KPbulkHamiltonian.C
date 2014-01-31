// $Id$

#include "KPbulkHamiltonian.h"
#include "Material.h"
#include "Messages.h"
#include "Constants.h"

using namespace Constants;



using namespace std;
//==================================================================


void KPbulkHamiltonian::do_print_info(void)
{

  ostringstream os;
  os << model_name<<"  ";
  os << "L = " << par.L1 << " M = " << par.M1 << " N = " << par.N1;

  if (model_name == "8x8")
  {
    os << "P0 = " << par.P0 << std::endl;
  }
  
  if (model_name == "14x14")
  {
    os << "P1 = " << par.P1_c << std::endl;
  }

  Messages::info(os.str());
  /*
  std::cout<<"6x6 parameters:"<<std::endl;
  std::cout<<"L = "<<result.L1+0.5  <<std::endl;
  std::cout<<"M = "<<result.M1+0.5  <<std::endl;
  std::cout<<"N = "<<result.N1/2.0  <<std::endl;
  std::cout<<"M-N = "<<result.M1+0.5-result.N1/2.0   <<std::endl;
  std::cout<<"M+N = "<<result.M1+0.5+result.N1/2.0  <<std::endl;
  std::cout<<"L-N = "<<result.L1+0.5-result.N1/2.0  <<std::endl;
  std::cout<<"L+2N = "<<result.L1+0.5+result.N1/2.0  <<std::endl;
  std::cout<<"N+ = "<<result.N1_xy  <<std::endl;
  std::cout<<"N- = "<<result.N1_yx  <<std::endl;
  std::cout<<"M-(N-) = "<<result.M1+0.5-result.N1_yx  <<std::endl;
  std::cout<<"M+(N-) = "<<result.M1+0.5+result.N1_yx  <<std::endl;
  std::cout<<"L-(N+) = "<<result.L1+0.5-result.N1_xy  <<std::endl;
  std::cout<<"L+2(N+) = "<<result.L1+0.5+result.N1_xy  <<std::endl;
  std::cout<<"---------------------------------------"<<std::endl;

  std::cout<<"8x8 parameters:"<<std::endl;
  std::cout<<"M-(N-) = "<<result.M1+0.5-result.N1_yx  <<std::endl;
  std::cout<<"M+(N-) = "<<result.M1+0.5+result.N1_yx  <<std::endl;
  std::cout<<"L-(N+) = "<<result.L1+0.5-result.N1_xy  <<std::endl;
  std::cout<<"L+2(N+) = "<<result.L1+0.5+result.N1_xy  <<std::endl;
  std::cout<<"---------------------------------------"<<std::endl;

  std::cout<<"14x14 parameters:"<<std::endl;
  std::cout<<"L = "<<result.L1+0.5  <<std::endl;
  std::cout<<"M = "<<result.M1+0.5  <<std::endl;
  std::cout<<"N = "<<result.N1/2.0  <<std::endl;
  std::cout<<"N+ = "<<result.N1_xy  <<std::endl;
  std::cout<<"N- = "<<result.N1_yx  <<std::endl;
  std::cout<<"M-N = "<<result.M1+0.5-result.N1/2.0   <<std::endl;
  std::cout<<"M+N = "<<result.M1+0.5+result.N1/2.0  <<std::endl;
  std::cout<<"L-N = "<<result.L1+0.5-result.N1/2.0  <<std::endl;
  std::cout<<"L+2N = "<<result.L1+0.5+result.N1/2.0  <<std::endl;
  std::cout<<"---------------------------------------"<<std::endl;

  */


}

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
 
  par.P0 = 0.0;
  par.E_c1 = 0.0;
  par.P1_c = 0.0;
  par.P2_c = 0.0;

}
//================================================================
KPbulkHamiltonian::~KPbulkHamiltonian()
{
}


//================================================================

KPbulkHamiltonian::KPbulkHamiltonian(const ModelOptions& options)
  : EFAbulkHamiltonian(options),
    model_name("6x6"),
    semiconductor(NULL),
    band_min(2),
    band_max(7)
{
}


void KPbulkHamiltonian::prepare_submodels(void)
{
  assert(semiconductor == NULL);

  ModelOptions opt =  get_options();
  opt.delete_all_submodels();
  semiconductor = Semiconductor::create(get_material(), opt);
  add_submodel("semiconductor", semiconductor);
}


//=================================================================

void KPbulkHamiltonian::do_init()
{

  EFAbulkHamiltonian::do_init();

  nullify_parameters();

  model_name = get_option("model", "6x6");

  model_name = get_option("kp_model", model_name);

  kpVVtermSymmetric = get_option("kpVVtermSymmetric", false);

  kpCVtermSymmetric = get_option("kpCVtermSymmetric", true);


  nullify_parameters();

  if (model_name == "2x2")
  {
    band_min = 0;
    band_max = 1;
    num_bands = 2;
  }
  else if (model_name == "6x6")
  {
    band_min = 2;
    band_max = 7;
    num_bands = 8;
  }
  else if (model_name == "8x8")
  {
    band_min = 0;
    band_max = 7;
    num_bands = 8;
  }
  else if (model_name == "14x14")
  {
    band_min = 0;
    band_max = 13;
    num_bands = 14;
  }
  else
  {
    Messages::error("Wrong kp model: "+model_name);
  }



  short i1 = 0;

  for (short i = band_min; i <= band_max; i++)
  {
      kp_bands.push_back(i);
      kp_bands_map.insert(make_pair (i,i1) );
      i1++;
  }
  
  //(re)set semiconductor model
  semiconductor->set_kp_model(model_name);

  //prepare k.p parameter
  semiconductor->calculate_kp_params(par);

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
  //create block-matrix 8x8 or 14x14
  const vector<KPbulkHamiltonian::MatrixElement>  temp(num_bands);


  Ham.resize(num_bands,temp);
  //--------------------------------------------------

  for (short i = 0; i < num_bands; i++)
  {
    for (short j = 0; j < num_bands; j++)
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
  }

  if (num_bands == 2)
  {
    // CB = 0, VB = 1
    Ham[0][0].quad[0][0] = 0.5*par.s2;
    Ham[0][0].quad[1][1] = 0.5*par.s2;
    Ham[0][0].quad[2][2] = 0.5*par.s1;

    Ham[1][1].quad[0][0] = par.L1 + 0.5;
    Ham[1][1].quad[1][1] = par.M1 + 0.5;
    Ham[1][1].quad[2][2] = par.M2 + 0.5;

    //Ham[1][1].quad[0][0] = par.M1 + 0.5;
    //Ham[1][1].quad[1][1] = par.L1 + 0.5;
    //Ham[1][1].quad[2][2] = par.M2 + 0.5;

    //Ham[1][1].quad[0][0] = par.M3 + 0.5;
    //Ham[1][1].quad[1][1] = par.M3 + 0.5;
    //Ham[1][1].quad[2][2] = par.L2 + 0.5;

    //cv part
    if (kpCVtermSymmetric)
    {
      // we assume an s-like hole orbital
      Ham[0][1].linear_left[0] = par.P2 * 0.5 *  Complex(0.0, 1.0);
      Ham[0][1].linear_left[1] = par.P2 * 0.5 *  Complex(0.0, 1.0);
      Ham[0][1].linear_left[2] = par.P1 * 0.5 *  Complex(0.0, 1.0);

      Ham[0][1].linear_right[0] = Ham[0][1].linear_left[0];
      Ham[0][1].linear_right[1] = Ham[0][1].linear_left[1];
      Ham[0][1].linear_right[2] = Ham[0][1].linear_left[2];

      for (short i = 0; i < 3; i++)
      {
        Ham[1][0].linear_left[i] =  conj(Ham[0][1].linear_left[i]);
        Ham[1][0].linear_right[i] = conj(Ham[0][1].linear_right[i]);
      }
    }
    else
    {
      Ham[0][1].linear_left[0] = par.P2 *  Complex(0, 1.0);
      Ham[0][1].linear_left[1] = par.P2 *  Complex(0, 1.0);
      Ham[0][1].linear_left[2] = par.P1 *  Complex(0, 1.0);

      for (short i = 0; i < 3; i++)
      {
        Ham[1][0].linear_right[i] = conj(Ham[0][1].linear_left[i]);
      }
    }




    Ham[0][0].constant += par.E_c;
    Ham[1][1].constant += par.E_v;
  }
  else // (num_bands >= 8)
  {
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

    for (short i = 0; i< 3; i++)
      for (short j = 0; j< 3; j++)
      {
        Ham[3][2].quad[i][j] = conj(Ham[2][3].quad[j][i]);
        Ham[4][3].quad[i][j] = conj(Ham[3][4].quad[j][i]);
        Ham[4][2].quad[i][j] = conj(Ham[2][4].quad[j][i]);
      }

    //-----------------------------------------------------------------------//
    // Crystal field splitting
    Ham[2][2].constant +=  par.d1 ;
    Ham[3][3].constant +=  par.d1 ;

    //=========================================================================//
    //---------Valence part without spin-orbit---------------------------------//
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

      Ham[0][2].linear_right[0] = Ham[0][2].linear_left[0];
      Ham[0][3].linear_right[1] = Ham[0][3].linear_left[1];
      Ham[0][4].linear_right[2] = Ham[0][4].linear_left[2];

      Ham[1][5].linear_left[0] = par.P2* 0.5 *  Complex(0.0, 1.0);
      Ham[1][6].linear_left[1] = par.P2* 0.5 *  Complex(0.0, 1.0);
      Ham[1][7].linear_left[2] = par.P1* 0.5 *  Complex(0.0, 1.0);

      Ham[1][5].linear_right[0] = Ham[1][5].linear_left[0];
      Ham[1][6].linear_right[1] = Ham[1][6].linear_left[1];
      Ham[1][7].linear_right[2] = Ham[1][7].linear_left[2];

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
      Ham[0][2].linear_left[0]= par.P2 *  Complex(0, 1.0);
      Ham[1][5].linear_left[0]= par.P2 *  Complex(0, 1.0);
      Ham[0][3].linear_left[1]= par.P2 *  Complex(0, 1.0);
      Ham[1][6].linear_left[1]= par.P2 *  Complex(0, 1.0);
      Ham[0][4].linear_left[2]= par.P1 *  Complex(0, 1.0);
      Ham[1][7].linear_left[2]= par.P1 *  Complex(0, 1.0);

      for (short i = 0; i< 2; i++)
        for (short j = 2; j< 8; j++)
          for (short  i1 = 0; i1<=2; i1++)
            Ham[j][i].linear_right[i1] = conj(Ham[i][j].linear_left[i1]);

    }
    //-----------------------------------------------------------------------//
    //absolute band edges
    //conduction band
    Ham[0][0].constant +=  par.E_c;
    Ham[1][1].constant +=  par.E_c;

    //valence band
    for (short i = 2; i < 8; i++)
      Ham[i][i].constant += par.E_v; 

    //-----------------------------------------------------------------------//
    //------------------spin-orbit interaction------------------------//
    Ham[2][3].constant   +=  Complex(0.0,-1.0)*par.d2;
    Ham[3][2].constant   +=  Complex(0.0, 1.0)*par.d2;

    Ham[2][7].constant += par.d3;
    Ham[7][2].constant += par.d3;

    Ham[5][4].constant += -par.d3;
    Ham[4][5].constant += -par.d3;

    Ham[6][4].constant +=  Complex(0.0,-1.0)*par.d3;
    Ham[4][6].constant +=  Complex(0.0,1.0)*par.d3;

    Ham[6][5].constant +=  Complex(0.0,-1.0)*par.d2;
    Ham[5][6].constant +=  Complex(0.0,1.0)*par.d2;

    Ham[7][3].constant +=  Complex(0.0,1.0)*par.d3;
    Ham[3][7].constant +=  Complex(0.0,-1.0)*par.d3;

    //-----------------------------------------------------------------!
    //==================================================================
    // 14x14 bands part
    if (num_bands > 8)
    {
      //second conduction band (Gamma7,8_c)
      for (short i = 8; i < 14; i++)
      {
        Ham[i][i].constant += par.E_c1;
        Ham[i][i].quad[0][0] = 0.5*par.s3;
        Ham[i][i].quad[1][1] = 0.5*par.s3;
        Ham[i][i].quad[2][2] = 0.5*par.s4;
      }

      // spin-orbit in second conduction

      Ham[8][9].constant   +=  Complex(0.0,-1.0)*par.d4;
      Ham[9][8].constant   +=  Complex(0.0, 1.0)*par.d4;

      Ham[8][13].constant += par.d4;
      Ham[13][8].constant += par.d4;

      Ham[11][10].constant += -par.d4;
      Ham[10][11].constant += -par.d4;

      Ham[12][10].constant +=  Complex(0.0,-1.0)*par.d4;
      Ham[10][12].constant +=  Complex(0.0,1.0)*par.d4;

      Ham[12][11].constant +=  Complex(0.0,-1.0)*par.d4;
      Ham[11][12].constant +=  Complex(0.0,1.0)*par.d4;

      Ham[13][9].constant +=  Complex(0.0,1.0)*par.d4;
      Ham[9][13].constant +=  Complex(0.0,-1.0)*par.d4;

      // spin-orbit in CV (cf in Tomic paper)

      Ham[3][8].constant   +=  Complex(0.0, -1.0)*par.d5;
      Ham[8][3].constant   +=  Complex(0.0, 1.0)*par.d5;

      Ham[2][9].constant   +=  Complex(0.0, 1.0)*par.d5;
      Ham[9][2].constant   +=  Complex(0.0, -1.0)*par.d5;

      Ham[6][11].constant   +=  Complex(0.0, 1.0)*par.d5;
      Ham[11][6].constant   +=  Complex(0.0, -1.0)*par.d5;

      Ham[5][12].constant   +=  Complex(0.0, -1.0)*par.d5;
      Ham[12][5].constant   +=  Complex(0.0, 1.0)*par.d5;

      Ham[8][7].constant   +=  -par.d5;
      Ham[7][8].constant   +=  -par.d5;

      Ham[9][7].constant   +=  Complex(0.0, 1.0)*par.d5;
      Ham[7][9].constant   +=  Complex(0.0, -1.0)*par.d5;

      Ham[10][5].constant   +=  par.d5;
      Ham[5][10].constant   +=  par.d5;

      Ham[10][6].constant   +=  Complex(0.0, -1.0)*par.d5;
      Ham[6][10].constant   +=  Complex(0.0, 1.0)*par.d5;

      Ham[2][13].constant   +=  -par.d5;
      Ham[13][2].constant   +=  -par.d5;

      Ham[3][13].constant   +=  Complex(0.0, 1.0)*par.d5;
      Ham[13][3].constant   +=  Complex(0.0, -1.0)*par.d5;

      Ham[4][11].constant   +=  par.d5;
      Ham[11][4].constant   +=  par.d5;

      Ham[4][12].constant   +=  Complex(0.0, -1.0)*par.d5;
      Ham[12][4].constant   +=  Complex(0.0, 1.0)*par.d5;

      // -----------------------------------------------------------------
      // CC couplings in symmetric form
      // up-up (_left)
      Ham[0][8 ].linear_left[0] = par.P1_c * 0.5 *  Complex(0.0, 1.0);
      Ham[0][9 ].linear_left[1] = par.P1_c * 0.5 *  Complex(0.0, 1.0);
      Ham[0][10].linear_left[2] = par.P1_c * 0.5 *  Complex(0.0,1.0);

      Ham[8][0].linear_left[0] = conj(Ham[0][8].linear_left[0]);
      Ham[9][0].linear_left[1] = conj(Ham[0][9].linear_left[1]);
      Ham[10][0].linear_left[2] = conj(Ham[0][10].linear_left[2]);

      // up-up (_right)
      Ham[0][8].linear_right[0] = Ham[0][8].linear_left[0];
      Ham[0][9].linear_right[1] = Ham[0][9].linear_left[1];
      Ham[0][10].linear_right[2] = Ham[0][10].linear_left[2];

      Ham[8][0].linear_right[0] = Ham[8][0].linear_left[0];
      Ham[9][0].linear_right[1] = Ham[9][0].linear_left[1];
      Ham[10][0].linear_right[2] = Ham[10][0].linear_left[2];

      // down-down (_left)
      Ham[1][11].linear_left[0] = par.P1_c * 0.5 *  Complex(0.0, 1.0);
      Ham[1][12].linear_left[1] = par.P1_c * 0.5 *  Complex(0.0, 1.0);
      Ham[1][13].linear_left[2] = par.P1_c * 0.5 *  Complex(0.0, 1.0);

      Ham[11][1].linear_left[0] = conj(Ham[1][11].linear_left[0]);
      Ham[12][1].linear_left[1] = conj(Ham[1][12].linear_left[1]);
      Ham[13][1].linear_left[2] = conj(Ham[1][13].linear_left[2]);

      // down-down (_right)
      Ham[11][1].linear_right[0] = Ham[11][1].linear_left[0];
      Ham[1][11].linear_right[0] = Ham[1][11].linear_left[0];

      Ham[12][1].linear_right[1] = Ham[12][1].linear_left[1];
      Ham[1][12].linear_right[1] = Ham[1][12].linear_left[1];

      Ham[13][1].linear_right[2] = Ham[13][1].linear_left[2];
      Ham[1][13].linear_right[2] = Ham[1][13].linear_left[2];

      //-----------------------------------------------------------------------
      // CV couplings in symmetric form
      // up-up (_left)
      if (kpCVtermSymmetric)
      {
        Ham[2][9 ].linear_left[2] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[2][10].linear_left[1] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[3][8 ].linear_left[2] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[3][10].linear_left[0] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[4][8 ].linear_left[1] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[4][9 ].linear_left[0] = par.P2_c * 0.5 *  Complex(0.0, -1.0);

        Ham[9][2].linear_left[2] = conj(Ham[2][9].linear_left[2]);
        Ham[10][2].linear_left[1] = conj(Ham[2][10].linear_left[1]);
        Ham[8][3].linear_left[2] = conj(Ham[3][8].linear_left[2]);
        Ham[10][3].linear_left[0] = conj(Ham[3][10].linear_left[0]);
        Ham[8][4].linear_left[1] = conj(Ham[4][8].linear_left[1]);
        Ham[9][4].linear_left[0] = conj(Ham[4][9].linear_left[0]);

        // up-up (_right)
        Ham[9][2].linear_right[2] = Ham[9][2].linear_left[2];
        Ham[2][9].linear_right[2] = Ham[2][9].linear_left[2];

        Ham[10][2].linear_right[1] = Ham[10][2].linear_left[1];
        Ham[2][10].linear_right[1] = Ham[2][10].linear_left[1];

        Ham[8][3].linear_right[2] = Ham[8][3].linear_left[2];
        Ham[3][8].linear_right[2] = Ham[3][8].linear_left[2];

        Ham[10][3].linear_right[0] = Ham[10][3].linear_left[0];
        Ham[3][10].linear_right[0] = Ham[3][10].linear_left[0];

        Ham[8][4].linear_right[1] = Ham[8][4].linear_left[1];
        Ham[4][8].linear_right[1] = Ham[4][8].linear_left[1];

        Ham[9][4].linear_right[0] = Ham[9][4].linear_left[0];
        Ham[4][9].linear_right[0] = Ham[4][9].linear_left[0];

        // down-down (_left)
        Ham[5][12].linear_left[2] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[5][13].linear_left[1] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[6][11].linear_left[2] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[6][13].linear_left[0] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[7][11].linear_left[1] = par.P2_c * 0.5 *  Complex(0.0, -1.0);
        Ham[7][12].linear_left[0] = par.P2_c * 0.5 *  Complex(0.0, -1.0);

        Ham[12][5].linear_left[2] = conj(Ham[5][12].linear_left[2]);
        Ham[13][5].linear_left[1] = conj(Ham[5][13].linear_left[1]);
        Ham[11][6].linear_left[2] = conj(Ham[6][11].linear_left[2]);
        Ham[13][6].linear_left[0] = conj(Ham[6][13].linear_left[0]);
        Ham[11][7].linear_left[1] = conj(Ham[7][11].linear_left[1]);
        Ham[12][7].linear_left[0] = conj(Ham[7][12].linear_left[0]);

        //down-down (_right)
        Ham[12][5].linear_right[2] = Ham[12][5].linear_left[2];
        Ham[5][12].linear_right[2] = Ham[5][12].linear_left[2];

        Ham[13][5].linear_right[1] = Ham[13][5].linear_left[1];
        Ham[5][13].linear_right[1] = Ham[5][13].linear_left[1];

        Ham[11][6].linear_right[2] = Ham[11][6].linear_left[2];
        Ham[6][11].linear_right[2] = Ham[6][11].linear_left[2];

        Ham[13][6].linear_right[0] = Ham[13][6].linear_left[0];
        Ham[6][13].linear_right[0] = Ham[6][13].linear_left[0];

        Ham[11][7].linear_right[1] = Ham[11][7].linear_left[1];
        Ham[7][11].linear_right[1] = Ham[7][11].linear_left[1];

        Ham[12][7].linear_right[0] = Ham[12][7].linear_left[0];
        Ham[7][12].linear_right[0] = Ham[7][12].linear_left[0];

        //--------------------------------------------------------
      }
      else
      {
        // up-up (_left)
        Ham[2][9 ].linear_left[2] = par.P2_c *  Complex(0.0, -1.0);
        Ham[2][10].linear_left[1] = par.P2_c *  Complex(0.0, -1.0);
        Ham[3][8 ].linear_left[2] = par.P2_c *  Complex(0.0, -1.0);
        Ham[3][10].linear_left[0] = par.P2_c *  Complex(0.0, -1.0);
        Ham[4][8 ].linear_left[1] = par.P2_c *  Complex(0.0, -1.0);
        Ham[4][9 ].linear_left[0] = par.P2_c *  Complex(0.0, -1.0);


        Ham[9 ][2].linear_right[2] = conj(Ham[2][9 ].linear_left[2]);
        Ham[10][2].linear_right[1] = conj(Ham[2][10].linear_left[1]);
        Ham[8 ][3].linear_right[2] = conj(Ham[3][8 ].linear_left[2]);
        Ham[10][3].linear_right[0] = conj(Ham[3][10].linear_left[0]);
        Ham[8 ][4].linear_right[1] = conj(Ham[4][8 ].linear_left[1]);
        Ham[9 ][4].linear_right[0] = conj(Ham[4][9 ].linear_left[0]);

        // down-down (_left)
        Ham[5][12].linear_left[2] = par.P2_c *  Complex(0.0, -1.0);
        Ham[5][13].linear_left[1] = par.P2_c *  Complex(0.0, -1.0);
        Ham[6][11].linear_left[2] = par.P2_c *  Complex(0.0, -1.0);
        Ham[6][13].linear_left[0] = par.P2_c *  Complex(0.0, -1.0);
        Ham[7][11].linear_left[1] = par.P2_c *  Complex(0.0, -1.0);
        Ham[7][12].linear_left[0] = par.P2_c *  Complex(0.0, -1.0);

        Ham[12][5].linear_right[2] = conj(Ham[5][12].linear_left[2]);
        Ham[13][5].linear_right[1] = conj(Ham[5][13].linear_left[1]);
        Ham[11][6].linear_right[2] = conj(Ham[6][11].linear_left[2]);
        Ham[13][6].linear_right[0] = conj(Ham[6][13].linear_left[0]);
        Ham[11][7].linear_right[1] = conj(Ham[7][11].linear_left[1]);
        Ham[12][7].linear_right[0] = conj(Ham[7][12].linear_left[0]);

      }
    }
  }


  //=====================================================================!
  //--------Valence band strain ------------------------------------
  //Ham[2][2].constant = par.l1s*strainM(1,1)+par.m1s*strainM(2,2)+par.m2s*strainM(3,3);
  //Ham[3][3].constant = par.m1s*strainM(1,1)+par.l1s*strainM(2,2)+par.m2s*strainM(3,3);
  //Ham[4][4].constant = par.m3s*strainM(1,1)+par.m3s*strainM(2,2)+par.l2s*strainM(3,3);
  //Ham[2][3].constant = par.n1s*strainM(2,1) ;
  //Ham[3][2].constant = par.n1s*strainM(2,1) ;
  //Ham[2][4].constant = par.n2s*strainM(3,1) ;
  //Ham[4][2].constant = par.n2s*strainM(3,1) ;
  //Ham[3][4].constant = par.n2s*strainM(3,2) ;
  //Ham[4][3].constant = par.n2s*strainM(3,2) ;


  //Transformation from crystal to calculation system
  for (short i = 0; i < num_bands; i++)
  {
    for (short j = 0; j < num_bands; j++)
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
  short nbands = num_bands;
  nbands = min(short(8), nbands);

  //nullify P matrix
  P.resize(3);
  for (short i = 0 ; i < 3; i++)
    {
      P[i].resize(nbands);
      for (short j = 0 ; j < nbands; j++)
        P[i][j].resize(nbands);
    }

  for (short pol = 0; pol < 3; pol++)
    for (short i = 0; i < nbands; i++)
      for (short j = 0; j < nbands; j++)
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





  //obtain P matrix (it is just on the 8x8 block)
  //we calculate derivative of H(k) matrix
  for (short band1 = 0; band1 < nbands; band1++)
    for (short band2 = 0; band2 < nbands; band2++)
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




void
KPbulkHamiltonian::get_hamiltonian_without_k(
    std::vector<std::vector<KPbulkHamiltonian::MatrixElement> >& ham) const
{
  ham.resize(band_max - band_min + 1);
  for (short i = 0; i <= band_max - band_min; i++)  ham[i].resize(band_max - band_min + 1);

  for (short i = 0; i <= band_max - band_min; i++)
    for (short j = 0; j <= band_max - band_min; j++)
      ham[i][j] = Ham[i + band_min][j + band_min];
}




//-------------------------------------------------------//

void KPbulkHamiltonian:: calculate_Hamiltonian_k_par (void)
{
  //allocation of the result
  //Initialization
   vector< vector<MatrixElement > > result = Ham;
  //--------------------------------------------------//


   for (short i = band_min; i <= band_max; i++)
   {
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
       {
         for (short j1 = 0; j1 < 3; j1++)
         {
           result[i][j].linear_left[i1]  += Ham[i][j].quad[i1][j1] * k_vector[j1];
           result[i][j].linear_right[j1] += Ham[i][j].quad[i1][j1] * k_vector[i1];
         }
       }

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



  const vector<double>  strain_Ham_Bir_Pikus1(num_bands,0.0);

  vector< vector<double> > strain_Ham_Bir_Pikus;
  strain_Ham_Bir_Pikus.resize(num_bands,strain_Ham_Bir_Pikus1);


  if (num_bands == 2)
  {
    //conduction band
    strain_Ham_Bir_Pikus[0][0]  = par.axs * ( strain_crystal(1,1) + strain_crystal(2,2) ) +
                                  par.azs * strain_crystal(3,3);

    strain_Ham_Bir_Pikus[1][1]  = par.l1s*strain_crystal(1,1)+
                                  par.m1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3);
    //strain_Ham_Bir_Pikus[1][1]  = par.m1s*strain_crystal(1,1)+
    //                              par.l1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3);
    //strain_Ham_Bir_Pikus[1][1]  = par.m3s*strain_crystal(1,1)+
    //                              par.m3s*strain_crystal(2,2)+par.l2s*strain_crystal(3,3);
  }
  else //if (num_bands > 2)
  {
  //---------------------------------------------------------------
  //  strain
  //---------------------------------------------------------------
  //valence band
  strain_Ham_Bir_Pikus[2][2]  = par.l1s*strain_crystal(1,1)+
                                par.m1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3);
  strain_Ham_Bir_Pikus[3][3]  = par.m1s*strain_crystal(1,1)+
                                par.l1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3);
  strain_Ham_Bir_Pikus[4][4]  = par.m3s*strain_crystal(1,1)+
                                par.m3s*strain_crystal(2,2)+par.l2s*strain_crystal(3,3);

  strain_Ham_Bir_Pikus[2][3]  = par.n1s*strain_crystal(2,1) ;
  strain_Ham_Bir_Pikus[3][2]  = par.n1s*strain_crystal(2,1) ;
  strain_Ham_Bir_Pikus[2][4]  = par.n2s*strain_crystal(3,1) ;
  strain_Ham_Bir_Pikus[4][2]  = par.n2s*strain_crystal(3,1) ;
  strain_Ham_Bir_Pikus[3][4]  = par.n2s*strain_crystal(3,2) ;
  strain_Ham_Bir_Pikus[4][3]  = par.n2s*strain_crystal(3,2) ;

  strain_Ham_Bir_Pikus[5][5]  = par.l1s*strain_crystal(1,1)+
                                par.m1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3);
  strain_Ham_Bir_Pikus[6][6]  = par.m1s*strain_crystal(1,1)+
                                par.l1s*strain_crystal(2,2)+par.m2s*strain_crystal(3,3);
  strain_Ham_Bir_Pikus[7][7]  = par.m3s*strain_crystal(1,1)+
                                par.m3s*strain_crystal(2,2)+par.l2s*strain_crystal(3,3);

  strain_Ham_Bir_Pikus[5][6]  = par.n1s*strain_crystal(2,1) ;
  strain_Ham_Bir_Pikus[6][5]  = par.n1s*strain_crystal(2,1) ;
  strain_Ham_Bir_Pikus[5][7]  = par.n2s*strain_crystal(3,1) ;
  strain_Ham_Bir_Pikus[7][5]  = par.n2s*strain_crystal(3,1) ;
  strain_Ham_Bir_Pikus[6][7]  = par.n2s*strain_crystal(3,2) ;
  strain_Ham_Bir_Pikus[7][6]  = par.n2s*strain_crystal(3,2) ;

  //conduction band
  strain_Ham_Bir_Pikus[0][0]  = par.axs * ( strain_crystal(1,1) + strain_crystal(2,2) ) +
                                par.azs * strain_crystal(3,3);
  strain_Ham_Bir_Pikus[1][1]  = par.axs * ( strain_crystal(1,1) + strain_crystal(2,2) ) + 
                                par.azs * strain_crystal(3,3);

  }


  //------------------------------------------------
  //potential
  //------------------------------------------------


  for (short i = 0; i < num_bands ; i++)
  {
    strain_Ham_Bir_Pikus[i][i] -= el_potential/Hartree;
  }


  //--------------------------------------------------
  //correction of the final Hamiltonian
  //--------------------------------------------------

  for (short i = 0; i <= band_max - band_min; i++)
    for (short j = 0; j <= band_max - band_min; j++)
      Hamiltonian[i][j].constant = Hamiltonian_without_strain_pot[i][j].constant +
                                   strain_Ham_Bir_Pikus[i + band_min][j + band_min];




}
//======================================================//

//===============================================================/
void KPbulkHamiltonian::calculate_optical_operator_k_par(void)
{
  short nbands = num_bands;
  nbands = min(short(8), nbands);

  vector <vector< vector<MatrixElement > > > P = P_gen;

   for (short i = 0; i < nbands; i++)
     for (short j = 0; j < nbands; j++)
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
