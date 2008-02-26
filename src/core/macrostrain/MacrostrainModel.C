#include "MacrostrainModel.h"
#include "Material.h"
#include <iostream>
MacrostrainModel::MacrostrainModel() : MacrostrainModelInterface() 
{
  stiffness = NULL;

  piezo     = NULL;

  poisson = NULL;

}


MacrostrainModel::~MacrostrainModel()
{

  PhysicalModelInterface::destroy(stiffness);
  PhysicalModelInterface::destroy(piezo);

}

//==========================================================================//

PhysicalModelInterface* MacrostrainModel::create_new (void) const
{
  return new MacrostrainModel();
}

//==========================================================================//
void MacrostrainModel::do_init()
{

  
  PhysicalModelInterface::destroy(stiffness);
  PhysicalModelInterface::destroy(piezo);

  const ModelOptions& opt =  get_options ();

  stiffness = Stiffness::create( get_material() -> get_structure(), opt  );
  
  stiffness->set_material(get_material());

  stiffness->init();


  piezo = Piezoelectricity::create( get_material() -> get_structure(), opt  ); 

  piezo->set_material(get_material());

  piezo->init();


  std::string poisson_name = opt.get_option("poisson_equation" , "no_poisson" );


  if (poisson_name != "no_poisson")
  {
    poisson =  SimulationInterface::find_simulation( poisson_name ) ;
    if (poisson == NULL)
      throw InitFailedException("MacrostrainModel: Unknown poisson model" + poisson_name);

    

    id_Ex = poisson->get_variable_id("Ex"); 
    Poisson_variables_ID.insert(id_Ex);
    
    id_Ey = poisson->get_variable_id("Ey"); 
    Poisson_variables_ID.insert(id_Ey);
    
    id_Ez = poisson->get_variable_id("Ez"); 
    Poisson_variables_ID.insert(id_Ez);

  
    

  }



}


//================================================================//
void MacrostrainModel::copy_from(const PhysicalModelInterface *rhs)
{
 
  // copy is not necessary as everything is created in do_init()
  
 
}

//==========================================================================//

void MacrostrainModel::read_database (void)
{
 
}

//==========================================================================//
void MacrostrainModel::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  const MacrostrainModel* matA = dynamic_cast< const MacrostrainModel*> (comp_A);

  const MacrostrainModel* matB = dynamic_cast< const MacrostrainModel*> (comp_B);

 

  stiffness->build_alloy(matA->stiffness, matB->stiffness, xa);
  
 

  piezo->build_alloy(matA->piezo, matB->piezo, xa);

  

  


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
    const std::set< ID > ids;
   
    std::vector< std::map< ID, double > >  field_components;

    bool got_solution = poisson->get_solution (element, q_point_vec, Poisson_variables_ID, field_components);

  

    if (got_solution)
    {
      Tensor1 field;
      

      field(1) = field_components[0][id_Ex];
      field(2) = field_components[0][id_Ey];
      field(3) = field_components[0][id_Ez];

     
      
      Material*   mat = get_material();

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


