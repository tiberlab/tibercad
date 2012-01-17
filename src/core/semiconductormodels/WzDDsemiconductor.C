// $Id$

#include "WzDDsemiconductor.h"
#include "WzSemiconductor.h"
#include "Constants.h"

using namespace std;
using namespace Constants;







//-------------------------------------------------------------/
void WzDDsemiconductor::do_calculate_conduction_band_extremum(void)
{
  vector<DDsemiconductor::band_extremum> result;

  const WzSemiconductor::WzDDparameters& par = (dynamic_cast<WzSemiconductor*> (semiconductor))->get_parameters ();


  double energy = par.Ev + par.EgGamma;

  if (strained) energy += (strain(1,1) + strain(2,2))* par.a_x + par.a_z *  strain(3,3);

  DDsemiconductor::band_extremum band_ext;
  band_ext.energy = energy;
  band_ext.degeneracy = 2;
  band_ext.mass_DOS = pow(par.m_c_xx * par.m_c_xx * par.m_c_zz, 1.0/3.0) ;

  result.push_back(band_ext);

  conduction_band = result;
}




