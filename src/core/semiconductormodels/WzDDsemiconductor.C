using namespace std;
#include "WzDDsemiconductor.h"
//--------------------------------------------------//
WzDDsemiconductor::WzDDsemiconductor(void): DDsemiconductor()
{
  EgGamma = 0;
  m_c_zz = 0; //conduction band mass along [0001] direction
  m_c_xx = 0; //conduction band mass along [1-210] direction
  A1 = 0;
  A2 = 0;
  A4 = 0;
  A5 =0 ;
  A6=0;
  A7=0; //k.p parameters
  a_x=0;
  a_z =0; //hydrostatic deformation potential
  D1 = 0;
  D2 = 0;
  D3 = 0;
  D4 = 0;
  D5 =0 ;
  D6 =0; //deformation potentials
  delta_s =0;
  delta_cr =0; //sp

}

//--------------------------------------------------//

WzDDsemiconductor::WzDDsemiconductor(const WzDDparameters& params): DDsemiconductor()
{
       Ev=params.Ev; //valence band energy
       EgGamma=params.EgGamma;//band gap
       m_c_zz=params.m_c_zz; //conduction band mass along [0001] direction
       m_c_xx=params.m_c_xx; //conduction band mass along [1-210] direction
       A1=params.A1;
       A2=params.A2;
       A3=params.A3;
       A4=params.A4;
       A5=params.A5;
       A6=params.A6;
       A7=params.A7; //k.p parameters
       a_x=params.a_x;
       a_z=params.a_z; //hydrostatic deformation potential
       D1=params.D1;
       D2=params.D2;
       D3=params.D3;
       D4=params.D4;
       D5=params.D5;
       D6=params.D6; //deformation potentials
       delta_s=params.delta_s;
       delta_cr=params.delta_cr; //spin
}
//---------------------------------------------------//
void WzDDsemiconductor::set_Eg(const double EgGamma1)
{
  EgGamma = EgGamma1;
}
//---------------------------------------------------//
void WzDDsemiconductor::set_mass_Gamma(const double m_c_zz_, const double m_c_xx_)
{
  m_c_zz = m_c_zz_;
  m_c_xx = m_c_xx_;
}
//----------------------------------------------------//
void WzDDsemiconductor::set_deform_pot_cond(const double a_x_, const double a_z_)
{
  a_x = a_x_;
  a_z = a_z_;
}
//----------------------------------------------------/
void WzDDsemiconductor::calculate_conduction_band_extremum(void)
{
  vector<DDsemiconductor::band_extremum> result;
  
  double energy = Ev + EgGamma;
 
  if (strained) energy += (strain(1,1) + strain(2,2))* a_x + a_z *  strain(3,3);
   
  DDsemiconductor::band_extremum band_ext;
  band_ext.energy = energy;
  band_ext.degeneracy = 2;
  band_ext.mass_DOS = pow(m_c_xx * m_c_xx * m_c_zz, 1/3) ;

  result.push_back(band_ext);

  conduction_band = result;
}
//----------------------------------------------------/

void  WzDDsemiconductor::calculate_valence_band_energy_extremum(void)
{


}
//----------------------------------------------------

