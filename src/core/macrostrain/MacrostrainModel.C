#include "MacrostrainModel.h"

MacrostrainModel::MacrostrainModel() : MacrostrainModelInterface() 
{
  stiffness = NULL;

  piezo     = NULL;

  crystal   = NULL;

}


MacrostrainModel::~MacrostrainModel()
{
  delete stiffness;

  delete piezo;

  delete  crystal;

}

//==========================================================================//

PhysicalModelInterface* MacrostrainModel::create_new (void) const
{
  return new MacrostrainModel();
}

//==========================================================================//
void MacrostrainModel::do_init()
{

  
  assert(stiffness != NULL);
  
  stiffness->init();


  assert(crystal != NULL);

  crystal->init();

  if (piezo != NULL) piezo->init();


}


//================================================================//
void MacrostrainModel::copy_from(const PhysicalModelInterface *rhs)
{
 
  const  MacrostrainModel*   temp = dynamic_cast<const MacrostrainModel* >  (rhs);

  
  stiffness = dynamic_cast<Stiffness* >( (temp->stiffness)->copy() );

  
  crystal= dynamic_cast<RotatedCrystal* >( (temp->crystal)->copy() );

  if (temp->piezo != NULL) piezo =  dynamic_cast<Piezoelectricity* >(  (temp->piezo)->copy() );
  

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
  crystal->build_alloy(matA->crystal, matB->crystal, xa);


  if (matA->piezo != NULL && matB->piezo != NULL)  piezo->build_alloy(matA->piezo, matB->piezo, xa);


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

//===========================================================================//
void MacrostrainModel::add_crystal(RotatedCrystal*  st)
{
  crystal = st;
}


//============================================================================//
