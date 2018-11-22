#ifndef _WZDDSEMICONDUCTOR_H_
#define _WZDDSEMICONDUCTOR_H_


#include "DDsemiconductor.h"
#include "PhysicalModel.h"

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

  

  //Destructor
  virtual ~WzDDsemiconductor(void) {};


  static WzDDsemiconductor* create(const ModelOptions& options);
 
 private:
 
  
  

 protected:

  //Constructor
  WzDDsemiconductor(const ModelOptions& options);

  PhysicalModel* create_new(void) const;
  
  //! calculates information about conduction bands
  /*!
    \f$ E_c^{\Gamma} = E_{c0}^{\Gamma} + a_{x}  (\varepsilon_{xx} + \varepsilon_{yy}) + a_z \varepsilon_{zz} \f$
  */
  virtual void  do_calculate_conduction_band_extremum(void);

};

inline
WzDDsemiconductor::WzDDsemiconductor(const ModelOptions& options)
 : DDsemiconductor(options)
{
}

inline PhysicalModel* WzDDsemiconductor::create_new( ) const
{
  return ( new WzDDsemiconductor(get_options()) );
}

inline WzDDsemiconductor* WzDDsemiconductor::create(const ModelOptions& options)
{
  return  new WzDDsemiconductor(options) ;
}





#endif 
