#ifndef ALLOY_H_
#define ALLOY_H_


#include <vector>
#include "Material.h"
#include "PhysicalProperties.h"
#include <point.h>

using  namespace  std;

class Alloy : public Material
{
 public:

  Alloy(const std::string& name, const std::string& structure = "zb") : Material(name, structure){};
	
    virtual ~Alloy();
	
    void init(const Dummy& database);
	
    void set_components(Material* matpoint);
	
	
  //  const  PhysicalProperties* Alloy::get_properties(const std::string& id,  const vector<double>& coord)  ; 
    const  PhysicalProperties*  get_properties(const std::string& id,  const Point& coord)  ; 
    
	
	
 private:

    vector<Material*>  components;  
    double  molar_fraction;
//	double calculate_molar_fraction(const vector<double>& coord);
	double calculate_molar_fraction(const Point& coord);
	
	
	
};

inline 
double 
Alloy::calculate_molar_fraction(const Point& coord) 
{
	return 0.0;
}




#endif /*ALLOY_H_*/
