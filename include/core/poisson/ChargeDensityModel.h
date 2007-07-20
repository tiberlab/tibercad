#ifndef _CHARGEDENSITYMODEL_H_
#define _CHARGEDENSITYMODEL_H_

#include "PhysicalModelInterface.h"
#include "SimulationInterface.h"
#include "elem.h"



//! This class computes the charge density
/*!


*/


class ChargeDensityModel : public PhysicalModelInterface
{

public:
  
  //!Constructor 
  ChargeDensityModel();
 
   //!Destructor
  ~ChargeDensityModel(){}

   //!provides electrons thermoelectric power [V/K]
  double  get_charge_density(void) const;
  
  static ChargeDensityModel* create();

 //!Recompute the charge density
  void re_init(void);

 //!Set the current element
  void set_element(const Elem* elem);

private:

  //! Current element 
  const Elem* _elem;

 
  //!Density charge
  double _charge_density; 
  
  //!A pointer to cherge density simulation
  SimulationInterface* _chd_sim;
  

protected:

  
  virtual void do_init (void);

  virtual void copy_from(const PhysicalModelInterface *rhs);

  virtual void read_database(void);

  virtual void read_bowing_parameters(void){};

  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 


  virtual PhysicalModelInterface* create_new (void) const;
 
 
 
};



inline
double   ChargeDensityModel::get_charge_density(void) const
{
  return _charge_density;
}


inline
ChargeDensityModel*  ChargeDensityModel::create()
{
  return (new  ChargeDensityModel());
}

inline
PhysicalModelInterface*  ChargeDensityModel::create_new () const
{
  return (new  ChargeDensityModel() ); 
}

inline
void 
ChargeDensityModel::set_element(const Elem* elem)
{

  _elem = elem;

}



#endif
