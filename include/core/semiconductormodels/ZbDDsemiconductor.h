// $Id$

#ifndef _ZBDDSEMICONDUCTOR_H_
#define _ZBDDSEMICONDUCTOR_H_
 
 
#include "DDsemiconductor.h"
#include "PhysicalModelInterface.h"
#include<vector>
#include<complex>

 
class ZbDDsemiconductor  : public DDsemiconductor
{
 public:
  
   
 
  
  virtual ~ZbDDsemiconductor(void) {};
 
  
 
  static ZbDDsemiconductor* create(const ModelOptions& options);
  
 private:
 
 
 
 protected:
 
  ZbDDsemiconductor(const ModelOptions& options) : DDsemiconductor(options) {};

  virtual PhysicalModelInterface* create_new(void) const;
   
  
  virtual void  do_calculate_conduction_band_extremum(void);
   
  
};
 
 
 
inline PhysicalModelInterface* ZbDDsemiconductor::create_new( ) const
{
  return ( new ZbDDsemiconductor(get_options()) );
}
 
inline ZbDDsemiconductor* ZbDDsemiconductor::create(const ModelOptions& options)
{
  return new ZbDDsemiconductor(options);
}
 
#endif


