#ifndef _WZDDSEMICONDUCTOR_H_
#define _WZDDSEMICONDUCTOR_H_


#include "DDsemiconductor.h"
#include "PhysicalModelInterface.h"

#include<vector>

//! A class to provide all neccessary parameters for drift-diffusion calculation for a wurtzite  crystal.
class  WzDDsemiconductor : public DDsemiconductor
/*!
  The class can calculate information about the band structure, such as
  band edge energy, effective mass for the density of states calculation and
  degeneracy.
  Conduction band masses do not depend on strain.
  Only \f$ \Gamma \f$ minimum of conduction band is considered. 
*/
{

 public:

  

  //Constructor
  WzDDsemiconductor(void) {};

  //Destructor
  ~WzDDsemiconductor(void) {};

 
  //! calculates information about conduction bands
  /*!
    \f$ E_c^{\Gamma} = E_{c0}^{\Gamma} + a_{x}  (\varepsilon_{xx} + \varepsilon_{yy}) + a_z \varepsilon_{zz} \f$
  */
  virtual void  calculate_conduction_band_extremum(void);

  //!calculates information about valence bands
  /*!
    Uses 6 band Luttinger kp theory 
  */
  virtual void  calculate_valence_band_extremum(void);

  

 

 

 
  static WzDDsemiconductor* create(void); 
 
 private:
 
  
  

 protected:
  PhysicalModelInterface* WzDDsemiconductor::create_new(void) const;
  


};

inline PhysicalModelInterface* WzDDsemiconductor::create_new( ) const
{
  return ( new WzDDsemiconductor() );
}

inline WzDDsemiconductor* WzDDsemiconductor::create() 
{
  return  new WzDDsemiconductor() ;
}





#endif 
