// $Id$


#include "DielectricModel.h"
#include "Material.h"
#include "Database.h"
#include "getpot.h"
//-------------------------------------------------------------------------//

DielectricModel::DielectricModel(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _dielectric_constant(0)
{
}


//-------------------------------------------------------------------------//


void  DielectricModel::do_init_alloy (const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa)
{
  const DielectricModel* modA = dynamic_cast<const DielectricModel*>(comp_A);

  const DielectricModel* modB = dynamic_cast<const DielectricModel*>(comp_B);

  alloy(_dielectric_constant,modA->_dielectric_constant, modB->_dielectric_constant,xa);


}

//-------------------------------------------------------------------------//

//--------------------------------------------------------//
void   DielectricModel::read_database(void)
{

  const Database& db = get_database();
  db.set_section("permittivity");
  _ep_x = db.get("permittivity", _ep_x);


}

//---------------------------------------------------------//

void  DielectricModel::do_init(void)
{


    ModelOptions& options = get_options ();

    if (has_option("permittivity"))
    {
      _ep_x = options.get_option("permittivity", _ep_x ) * 1e-2 ;

      _dielectric_constant(1,1) = _ep_x;

      _dielectric_constant(2,2)= _ep_x;

      _dielectric_constant(3,3)= _ep_x;
    }

}






//-------------------------------------------------------------------------//
