#include "MacrostrainModel.h"
#include "Material.h"
MacrostrainModel::MacrostrainModel() : MacrostrainModelInterface() 
{
  stiffness = NULL;

  piezo     = NULL;

  

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

  
 

  const ModelOptions& opt =  get_options ();

  stiffness = Stiffness::create( get_material() -> get_structure(), opt  );
  
  stiffness->set_material(get_material());

  stiffness->init();


  piezo = Piezoelectricity::create( get_material() -> get_structure(), opt  ); 

  piezo->set_material(get_material());

  piezo->init();


}


//================================================================//
void MacrostrainModel::copy_from(const PhysicalModelInterface *rhs)
{
 
  // copy is not necessary as they are created in do_init()
  
  //const  MacrostrainModel*   temp = dynamic_cast<const MacrostrainModel* >  (rhs);
  
  //stiffness = dynamic_cast<Stiffness* >( (temp->stiffness)->copy() );

  //piezo =  dynamic_cast<Piezoelectricity* >(  (temp->piezo)->copy() );
  

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
