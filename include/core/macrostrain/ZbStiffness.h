#ifndef _ZBSTIFFNESS_H_
#define _ZBSTIFFNESS_H_

#include "Stiffness.h"

/*! assembles stiffness tensor in crystal system for a zinc-blende crystal
  \f$ c_{11} = C_{xxxx}, c_{12}=C_{xxyy}, c_{44}=C_{xyxy} \f$ (Voigt notation)
*/
class ZbStiffness : public Stiffness
{

 public:

  //!method that sets Young moduli
  void set_moduli(double c11, double c12, double c44);


  static ZbStiffness* create(const ModelOptions& options);

 protected:

  //!Empty constructor
  ZbStiffness(const ModelOptions& options);

  //reads database
  void read_database ( );

  virtual PhysicalModel* create_new(void) const;

  virtual void do_init (void);


};


inline  PhysicalModel* ZbStiffness::create_new(void) const
{
  return ( new ZbStiffness(get_options()) ) ;
}


inline ZbStiffness* ZbStiffness::create(const ModelOptions& options)
{
   return new ZbStiffness(options) ;
}

#endif
