// $Id$

#include "MacrostrainModel.h"
#include "Material.h"
#include <iostream>



MacrostrainModel::MacrostrainModel(const ModelOptions& options)
 : MacrostrainModelInterface(options)
{
  stiffness = NULL;

  piezo     = NULL;

  poisson = NULL;

}


MacrostrainModel::~MacrostrainModel()
{
}

//==========================================================================//

PhysicalModelInterface* MacrostrainModel::create_new (void) const
{
  return new MacrostrainModel(get_options());
}



void
MacrostrainModel::prepare_submodels(void)
{
  assert(stiffness == NULL);
  assert(piezo == NULL);

  ModelOptions opt(get_options());
  opt.delete_all_submodels();
  stiffness = Stiffness::create(opt);
  add_submodel("stiffness", stiffness);

  piezo = Piezoelectricity::create(get_material(), opt);
  add_submodel("piezo", piezo);
}


//==========================================================================//
void MacrostrainModel::do_init()
{

  const ModelOptions& opt =  get_options ();

  std::string poisson_name = opt.get_option("poisson_equation" , "no_poisson" );

  if (poisson_name != "no_poisson")
  {
    poisson =  SimulationInterface::find_simulation( poisson_name ) ;
    if (poisson == NULL)
      throw InitFailedException("MacrostrainModel: Unknown poisson model" + poisson_name);

    id_E = poisson->get_solution_id("ElField");
  }

}



//==========================================================================//
void MacrostrainModel::add_stiffness(Stiffness*  st)
{
  stiffness = st;
}

//==========================================================================//
void MacrostrainModel::add_piezo(Piezoelectricity* pz)
{
  piezo = pz;
}



//============================================================================//
void MacrostrainModel::get_converse_piezo_stress(Tensor2Sym& sigma, const Elem* element)
{

  //crystal system
  sigma = 0;

  if (poisson != NULL)
  {

    Point q_point = element->centroid();
    std::vector<Point> q_point_vec(1, q_point);

    std::map<ID, std::vector<double> > field_components;
    field_components[id_E].resize(0);

    bool got_solution = poisson->get_solution(element, field_components, q_point_vec);

    if (got_solution)
    {
      Tensor1 field;

      field(1) = field_components[id_E][0];
      field(2) = field_components[id_E][1];
      field(3) = field_components[id_E][2];

      const Material*   mat = get_material();

      const RotatedCrystal&   cr = mat->get_rotated_crystal ();

      field = cr.RotMatrix.transpose() * field; //convert to crystal system [V/cm];

      field =  100 * field ; //convert [V/m]

      field = field / Constants::field_gauss_unit; //convert to gauss units


      piezo->calculate_product_by_vector(field, sigma);//calculate in the crystal system




      sigma =  sym( cr.RotMatrix * sigma  * cr.RotMatrix.transpose()); //convert to calculation system



      // sigma = sigma / Constants::polarization_gauss_unit; //SI units

      sigma = sigma * Constants::c * 10.0 / 1e4 / 10.0; //SI units

      //  std::cerr << Constants::polarization_gauss_unit << "\n";

      sigma = sigma * 1e-9; //GPa units


    }
  }

}


