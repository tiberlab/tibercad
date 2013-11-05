// $Id$

#include "SBZbCondBandBulkHamiltonian.h"
#include "ZbSemiconductor.h"
#include "Constants.h"
#include "Messages.h"

using namespace std;
using namespace Constants;


void 
SBZbCondBandBulkHamiltonian::do_print_info(void)
{

  ostringstream os;
  os << "single-band  ";

  os << "Mass = " << zb_par->m_G << std::endl;

  //os << "M_xx= " << imass(1,1)  << "M_yy= " << imass(2,2) 
  //   << "M_zz= " << imass(3,3) << std::endl;
  //os << "M_xy= " << imass(1,2)  << "M_xz= " << imass(1,3) 
  //   << "M_yz= " << imass(2,3) << std::endl;

  Messages::info(os.str());
}



//===========================================================================//
void SBZbCondBandBulkHamiltonian::calculate_for_init(void)
{
  const ZbSemiconductor::ZbDDparameters& par = (dynamic_cast<ZbSemiconductor*> (semiconductor)) -> get_parameters();

  zb_par = &par;

  //std::cout<<"calculate_for_init"<<std::endl;

  const ModelOptions& opt =  get_options ();

  min_name = opt.get_option("minimum_name", "Gamma");


  if (min_name == "Gamma")
    {
      imass = Tensor2Sym(0);
      imass(1,1) = 1.0/par.m_G;
      imass(2,2) = 1.0/par.m_G;
      imass(3,3) = 1.0/par.m_G;

      edge = (par.Ev  + par.EgGamma) / Hartree;

      kp_bands.resize(1,0);

    }

  if ( (min_name == "X") || (min_name == "x") )
    {
      if ( opt.get_option("valley", 1) == 1)
	{
	  imass = Tensor2Sym(0);
	  imass(1,1) = 1.0/par.m_l_X;
	  imass(2,2) = 1.0/par.m_t_X;
	  imass(3,3) = 1.0/par.m_t_X;

	  edge = (par.Ev  + par.EgX) / Hartree;

	  kp_bands.resize(1,0);
	}
      if ( opt.get_option("valley", 1) == 2)
	{
	  imass = Tensor2Sym(0);
	  imass(1,1) = 1.0/par.m_t_X;
	  imass(2,2) = 1.0/par.m_l_X;
	  imass(3,3) = 1.0/par.m_t_X;

	  edge = (par.Ev  + par.EgX) / Hartree;

	  kp_bands.resize(1,0);
	}
      if ( opt.get_option("valley", 1) == 3)
	{
	  imass = Tensor2Sym(0);
	  imass(1,1) = 1.0/par.m_t_X;
	  imass(2,2) = 1.0/par.m_t_X;
	  imass(3,3) = 1.0/par.m_l_X;

	  edge = (par.Ev  + par.EgX) / Hartree;

	  kp_bands.resize(1,0);
	}

        _degeneracy = 2;
    }

  calculate_Hamiltonian_gen();

  calculate_Hamiltonian_k_par();



}

//===========================================================================//
void SBZbCondBandBulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{

  //TODO: Redundance with methods in ZbDDSemiconductor::do_calculate_conduction_band_extremum.
  //Put anything in only one method
  const ZbSemiconductor::ZbDDparameters& par = (dynamic_cast<ZbSemiconductor*> (semiconductor)) -> get_parameters();

  zb_par = &par;

  //now strain and potential
  if (min_name == "Gamma")
    {
       Hamiltonian[0][0].constant =  Hamiltonian_without_strain_pot[0][0].constant - el_potential/Hartree
	 + (zb_par->a_c) * trace(strain_crystal )/Hartree;

    }

  else if (min_name == "X")
  {
    //Hydrostatic strain
    double hydro_def_pot = 0.0;
    double uniaxial_def_pot = 0.0;
    double def_pot = 0.0;

    hydro_def_pot = (par.def_vol_X - (1.0/3.0)*par.def_uniax_X) * trace(strain_crystal);

    //Uniaxial strain part----------------------------------
          Tensor1 k;
          vector<double> dE_uniax(3);

          k(1) = 1 ; k(2) =  0; k(3) =  0;
          dE_uniax[0] = par.def_uniax_X * (k * (strain_crystal * k));

          k(1) = 0 ; k(2) =  1; k(3) =  0;
          dE_uniax[1] = par.def_uniax_X * (k * (strain_crystal * k));

          k(1) = 0 ; k(2) =  0; k(3) =  1;
          dE_uniax[2] = par.def_uniax_X * (k * (strain_crystal * k));

          for (short i = 0; i <3; i++)
          {
            uniaxial_def_pot += dE_uniax[i];
          }

          def_pot = uniaxial_def_pot + hydro_def_pot;

          Hamiltonian[0][0].constant =  Hamiltonian_without_strain_pot[0][0].constant - el_potential/Hartree
                   + def_pot/Hartree;

  }

}
