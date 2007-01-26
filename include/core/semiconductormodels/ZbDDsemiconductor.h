#ifndef _ZBDDSEMICONDUCTOR_H_
#define _ZBDDSEMICONDUCTOR_H_


#include "DDsemiconductor.h"
#include "PhysicalModelInterface.h"
#include<vector>
#include<complex>


//! A class to provide all neccessary parameters for drift-diffusion calculation for a zinc-blend (or diamond) crystal
class ZbDDsemiconductor  : public DDsemiconductor
/*!
  The class can calculate information about the band structure, such as
  band edge energy, effective mass for the density of states calculation and
  degeneracy.
  Conduction band masses do not depend on strain.
*/
{
 public:
 
  


  //!Constructor  
  ZbDDsemiconductor(void) {};

  //!Destructor
  ~ZbDDsemiconductor(void) {};

  //! calculates information about conduction bands
  /*!
    \f$ E_c^{\Gamma} = E_{c0}^{\Gamma} + a_c  \mathop{\rm Tr} (\varepsilon_{ij}) ;\f$

    \f$ E_c^{X} = E_{c0}^{X} + \Xi_v^{X}  \mathop{\rm Tr} (\varepsilon_{ij}) +
     \Xi_u^{X} k_i k_j \varepsilon_{ij} ;\f$

    \f$ E_c^{L} = E_{c0}^{L} + \Xi_v^{L}  \mathop{\rm Tr} (\varepsilon_{ij}) +
     \Xi_u^{L} k_i k_j \varepsilon_{ij} ,\f$

    where k is a unit vector from the Brillouin zone center to a minimum point 
  */
  virtual void  calculate_conduction_band_extremum(void);
  
  //! calculates information about valence bands
  /*!
    Uses 6 band Luttinger kp theory 
  */
  virtual void  calculate_valence_band_extremum(void);

  static ZbDDsemiconductor* create();

 private:



 protected:

  virtual PhysicalModelInterface* create_new(void) const;

 
  

 
};



inline PhysicalModelInterface* ZbDDsemiconductor::create_new( ) const
{
  return ( new ZbDDsemiconductor() );
}

inline ZbDDsemiconductor* ZbDDsemiconductor::create()
{
  return new ZbDDsemiconductor();
}

#endif 
