
//

//
//
#ifndef EMPROPERTIES_H
#define EMPROPERTIES_H

#include "PhysicalProperties.h"
#include "Material.h"
//#include "DataBaseCall.h"
// forward declarations
class Dummy;

class EmProperties : public PhysicalProperties
{
public:
  // Constructor calls base class  constructor with  parameter 'id'  and  then empty constructor of  derived class is  executed
  EmProperties(const std::string& id):PhysicalProperties( id ) {};


  ~EmProperties() {};

  virtual void read_database_bowing_parameters(const Dummy& db);
  //  virtual void read_database(DataBaseCall& db);
 virtual void read_database(Dummy& db);


  void set_properties_alloy(const  PhysicalProperties* prop_comp1,
                            const  PhysicalProperties* prop_comp2,
                            double molar_fraction) ;

  double  get_lattice_constant_a();

private:
  double latt_const_a;


};



#endif


