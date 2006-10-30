#include "Alloy.h"

//Alloy::Alloy(): Material(const std::string& name, const std::string& structure = "zb")
//{
//}

Alloy::~Alloy()
{
}


void 
Alloy::set_components(Material* matpoint)
{	
  //  vector<Material*> components .push_back(matpoint)
  components.push_back(matpoint);
		
}


//virtual :  for  ALLOY  add properties object also for components.
void 
Alloy::add_properties(PhysicalProperties* properties)

{
  const std::string& id = properties->get_id();
  _properties[id] = properties;
  
  cout << " in  Alloy::add_properties" << endl;
  
  //components[0]->add_properties(properties);
  cout << "components[0]->get_name()  "   << components[0]->get_name() << endl;
  
  
  
  
  //components[1]->add_properties(properties);
  cout << "components[1]->get_name()  "   << components[1]->get_name() << endl;
  
}



// for Alloy  object: virtual init writes in PhysicalProperties only the  bowing  parameters
//  AND  call init  for  component materials  of  the  alloy.
//
void
//Alloy::init(DataBaseCall& database)
Alloy::init(const Dummy& database)

{
  PropertyMap::iterator it = _properties.begin();
  const PropertyMap::const_iterator end = _properties.end();

  for ( ; it != end; ++it)
  {
    (it->second)->set_material(this);
   // (it->second)->read_database_bowing_parameters(database);
  }
    
  //components[0]->init(database);
    
  cout << "components[0]->_name  " << components[0]->get_name() << endl;
    
  //components[1]->init(database);
  cout << "components[1]->_name  " << components[1]->get_name()<<  endl;
    
    
}






// virtual  in base class (Material)
const  PhysicalProperties* 
Alloy::get_properties(const std::string& id,  const Point& coord) 
{		
  // calculates a possibly position_dependent molar fraction x =f(r)
  molar_fraction = calculate_molar_fraction(coord);
  cout << " molar_fraction  in get_properties******* = " << molar_fraction << endl;
   
   
  const PropertyMap::const_iterator end = _properties.end();
  PropertyMap::const_iterator it = _properties.find(id);
  if (it != end)
  { // do for the PhysicalProperties it->second
   	  	
    const PhysicalProperties* prop_comp1 = components[0]->get_properties(id, coord);
    cout << "prop_comp1 ->latt_const_a" << ((EmProperties*)prop_comp1) ->get_lattice_constant_a()<< endl;
      
      
    // get properties object  from alloy comp 1
    const PhysicalProperties* prop_comp2 = components[1]->get_properties(id, coord);
    cout << "prop_comp2 ->latt_const_a" << ((EmProperties*)prop_comp2) ->get_lattice_constant_a() << endl;
    
    (it->second) -> set_properties_alloy(prop_comp1, prop_comp2, molar_fraction);
    return it->second;
  }
  else
    return NULL;
      	
}
	
	
	




















