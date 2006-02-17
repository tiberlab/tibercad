#include "Alloy.h"

//Alloy::Alloy(): Material(const std::string& name, const std::string& structure = "zb")
//{
//}

Alloy::~Alloy()
{
}


void Alloy::set_components(Material* matpoint)

{
	
  //  vector<Material*> components .push_back(matpoint)
  components.push_back(matpoint);
	
	
}


// for Alloy  object: virtual init writes in PhysicalProperties only the  bowing  parameters
void
Alloy::init(const Dummy& database)
{
  PropertyMap::iterator it = _properties.begin();
  const PropertyMap::const_iterator end = _properties.end();

  for ( ; it != end; ++it)
    {
      (it->second)->set_material(this);
      (it->second)->read_database_bowing_parameters(database);
    }
}



// virtual  in base class (Material)
const  PhysicalProperties* Alloy::get_properties(const std::string& id,  const Point& coord) 

{
	
	
  // calculates a possibly position_dependent molar fraction x =f(r)
  molar_fraction = calculate_molar_fraction(coord);
	
  const PropertyMap::const_iterator end = _properties.end();
  PropertyMap::const_iterator it = _properties.find(id);
  if (it != end)
    { // do for the PhysicalProperties it->second
   	  	
      const PhysicalProperties* prop_comp1 = components[1]->get_properties(id, coord);
      // get properties object  from alloy comp 1
      const PhysicalProperties* prop_comp2 = components[2]->get_properties(id, coord);
    
      (it->second) -> set_properties_alloy(prop_comp1, prop_comp2, molar_fraction);
      return it->second;
    }
  else
    return NULL;
    
  
	
}
	
	
	




















