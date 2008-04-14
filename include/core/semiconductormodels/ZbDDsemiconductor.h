#ifndef _ZBDDSEMICONDUCTOR_H_
#define _ZBDDSEMICONDUCTOR_H_
 
 
#include "DDsemiconductor.h"
#include "PhysicalModelInterface.h"
#include<vector>
#include<complex>

 
class ZbDDsemiconductor  : public DDsemiconductor
{
 public:
  
   
 
 
  ZbDDsemiconductor(void) {};
  
  virtual ~ZbDDsemiconductor(void) {};
 
  
 
  static ZbDDsemiconductor* create();
  
 private:
 
 
 
 protected:
 
  virtual PhysicalModelInterface* create_new(void) const;
   
  
 
  virtual void  do_calculate_conduction_band_extremum(void);
   
  
  virtual void  do_calculate_valence_band_extremum(void);
 
  
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


