#ifndef _WZSTIFFNESS_H_
#define _WZSTIFFNESS_H_

#include "Stiffness.h"
#include "PhysicalModelInterface.h"
class WzStiffness : public Stiffness
{

 public:

  //!Empty constructor
  WzStiffness();
  
  //!Constructor that sets moduli

  WzStiffness(double c11, double c12, double c13, double c33, double c44 );


  //! method that sets moduli
  /*!assembles stiffness tensor in crystal system for a wurtzite crystal
    following Eq (16) J. of Physics Condensed Matter v 14 p.3399 O. Ambacher et al 
  z axis is parallel to [0001] direction
  (x,y,z) - othogonormal system 
   \f$ c_{11} = C_{xxxx}, c_{12}=C_{xxyy}, c_{13} = C_{xxzz}, c_{33} = C_{zzzz}, c_{44} = C_{yzxz} \f$
  */  
  void set_moduli(double c11, double c12, double c13, double c33, double c44 );
  

  static WzStiffness* create(void);

 protected:
  
 
  virtual void read_database( );


  virtual PhysicalModelInterface* create_new(void) const;

  virtual void do_init(void);
  
};


inline PhysicalModelInterface* WzStiffness::create_new(void) const
{
  return ( new WzStiffness() );
}


inline WzStiffness* WzStiffness::create()
{
  return new WzStiffness();
}
#endif
