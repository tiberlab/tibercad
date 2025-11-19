/*
 * incidentwave.C
 *
 *  Created on: 4 Oct 2021
 *      Author: pamiri
 */

#include "Dipole_source.h"
#include "tibercad/module/TiberModule.h"

void Dipole_source::do_init(void)
{
  write_type("Dipole Source");
  get_parameter("kr_ratio", _kr);
  get_parameter("steps", _steps);
  set_dipole_elements(_kr,_steps);
 // typer = "Incident Wave";
}
void Dipole_source::Calculate_M_Matrix(void)
{
  ;;
}

