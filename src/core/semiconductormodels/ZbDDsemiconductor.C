// $Id$

#include "ZbDDsemiconductor.h"
#include "ZbSemiconductor.h"
#include "Semiconductor.h"
#include "Constants.h"

using namespace std;
using namespace Constants;


 //---------------------------------------------//


 void ZbDDsemiconductor::do_calculate_conduction_band_extremum(void)

 {
    vector<DDsemiconductor::band_extremum>   result;

    const ZbSemiconductor::ZbDDparameters& par = (dynamic_cast<ZbSemiconductor*> (semiconductor))->get_parameters ();

    band_extremum  band_ex;

    vector<string> valleys = {"Gamma", "X", "L"};
    get_option("valley", valleys);

    //
    // there are 3 types of extrema: Gamma, X and L
    //

    //---------------------------------------------------------
    // Gamma minimum
    //---------------------------------------------------------
    auto it = std::find(valleys.begin(), valleys.end(), "Gamma");
    if (it != valleys.end())
    {
      double Ec_G = par.Ev + par.EgGamma;

      if (strained)  Ec_G  += par.a_c * trace(strain);

      band_ex.degeneracy = 2   ; //spin
      band_ex.energy     = Ec_G;
      band_ex.mass_DOS   = par.m_G ; //mass is isotropic

      result.push_back(band_ex);
    }


    //----------------------------------------------------------------
    // L minimum
    //----------------------------------------------------------------
    it = std::find(valleys.begin(), valleys.end(), "L");
    if (it != valleys.end())
    {
      double Ec_L = par.Ev + par.EgL;

      if (strained)
      {

        //Hydrostatic strain part--------------------------------
        Ec_L  += (par.def_vol_L - (1.0/3.0)*par.def_uniax_L ) * trace(strain);
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
        //band_ex.degeneracy = 16   ;
        /* spin degeneracy and  8 equivalent minima (
           [1,1, 1],[-1,-1,-1],
           [-1,1,1],[1,-1,-1]
           [1,-1,1],[-1,1,-1]
           [1,1,-1],[-1,-1,1]
         */

        // Matthias 2015-09-24: let's always put them explicitly
        band_ex.degeneracy = 4   ;
        band_ex.energy     = Ec_L;
        band_ex.mass_DOS = pow(par.m_t_L * par.m_t_L *  par.m_l_L, 1.0/3.0 );
        for (short i = 0; i <4; i++)
          result.push_back(band_ex);

      }
    }

    //------------------------------------------------------------------
    // X minimum
    //------------------------------------------------------------------
    it = std::find(valleys.begin(), valleys.end(), "X");
    if (it != valleys.end())
    {
      double Ec_X = par.Ev + par.EgX;

      if (strained)
      {
        //Hydrostatic strain part--------------------------------
        Ec_X  += (par.def_vol_X - (1.0/3.0)*par.def_uniax_X) * trace(strain);
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
          band_ex.degeneracy = 4; //spin and k->-k
          band_ex.mass_DOS = mass_DOS;
          band_ex.energy = Ec_X + dE_uniax[i];
          result.push_back(band_ex);
        }


      }
      else
      {
        band_ex.energy     = Ec_X;
        //band_ex.degeneracy = 12   ; // spin degeneracy and 6 equivalent minima
        // Matthias 2015-09-24: let's always put them explicitly
        band_ex.degeneracy = 4;
        band_ex.mass_DOS = pow(par.m_t_X * par.m_t_X *  par.m_l_X, 1.0/3.0 );
        for (short i = 0; i < 3; i++)
          result.push_back(band_ex);
      }
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





